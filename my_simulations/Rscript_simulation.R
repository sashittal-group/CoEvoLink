# Script for simulation study performed in Braga et al. 2020 
# Bayesian inference of ancestral host-parasite interactions under a phylogenetic model of host repertoire evolution

library(ape)
library(tidyverse)
library(igraph)
library(MCMCpack)
require(phylotate)
require(data.table)
library(coda)
library(kdensity)


# read phylos
path_sim_data <- "./data/"

host_tree <- read.tree(paste0(path_sim_data,"angio_25tips_bl1.phy"))
hosts <- host_tree$tip.label

tree <- read.tree(paste0(path_sim_data,"Nymphalini.phy"))
phylob <- tree$tip.label

# Read all simulated datasets ----

data_sims <- tibble()

clocks <- c(0.1,0.5,1)     # values of overall rate (clock)
betas <- c(0,1,4)          # values of beta
iterations <- 50           # number of iterations/replicates

for(b in betas){
  for(c in clocks){
    for(i in 1:iterations){    
      name <- paste0("sim",i,"-b",b,"-c",c)
      data <- as.data.frame(read.nexus.data(paste0(path_sim_data,"sim",i,"-b",b,"-c",c,".nex")))
      #assign(name, data)
      rownames(data) <- hosts
        
      graph <- graph_from_incidence_matrix(t(data), weighted = TRUE)
      el <- get.data.frame(graph, what = "edges") %>% 
        mutate(beta = b, clock = c, it = i)
      
      data_sims <- rbind(data_sims, el)
    }
  }
}

internal_sims <- filter(data_sims, to %in% paste0("Index_",35:67)) %>% 
  complete(from, to, beta, clock, it, fill = list(weight = 0))
tips_sims <- filter(data_sims, to %in% phylob)


# Plot simulated data sets for one iteration of parameter combinations ----

plot_sims <- filter(tips_sims, it == 1) %>%
  mutate(to = factor(to, levels = phylob),
  from = factor(from, levels = hosts_ladder))

plot_sims %>% 
  mutate(beta_lab = factor(beta, labels = c(beta ~"= 0", beta ~"= 1", beta ~"= 4")),
         clock_lab = factor(clock, labels = c(mu ~"= 0.1", mu ~"= 0.5", mu ~"= 1.0"))) %>% 
  ggplot(aes(x = from, y = to, fill = weight)) +
  geom_tile() +
  scale_fill_manual(values = c("grey50", "black")) +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  facet_grid(clock_lab ~ beta_lab, labeller = label_parsed) +
  theme_bw() +
  labs(y = "Rate of repertoire evolution", x = "Effect of host phylogeny") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey40"),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    axis.ticks = element_blank(),
    strip.text.y = element_text(angle = 0),
    legend.position = "none")



# Inference ---------------------------------------------------
# Read log files

path_sim_out <- "./output/"

clocks <- c(0.1,0.5,1)     # values of overall rate (clock)
betas <- c(0,1,4)          # values of beta
iterations <- 50           # number of iterations/replicates

log_sim <- tibble()

for(b in betas){
  for(c in clocks){
    for(i in 1:iterations){
      for(chain in 1:2){                       # change here #
        post <- read.table(paste0(path_sim_out,"out.bl1-rep",chain,"-sim",i,"-b",b,"-c",c,".log"), header = TRUE)[,c(1,2,5,6,8:11)]
        colnames(post) <- c("generation","posterior","clock","beta", "0_to_1", "1_to_0", "1_to_2", "2_to_1")
        post <- mutate(post, true_beta = b, true_clock = c, iteration = i, chain = chain) 
        log_sim <- rbind(log_sim, post)
      }
    }
  }
}


# Convergence test ----

gelman_sim_thin <- tibble()

for(b in betas){
  for(c in clocks){
    for(i in 1:iterations){
      chain1 <- filter(log_sim, true_beta == b, true_clock == c, iteration == i, chain == 1) %>% 
        filter(generation %in% seq(10000,60000,200)) %>% 
        dplyr::select(clock, beta, contains("_to_"))
      chain2 <- filter(log_sim, true_beta == b, true_clock == c, iteration == i, chain == 2) %>%
        filter(generation %in% seq(10000,60000,200)) %>% 
        dplyr::select(clock, beta, contains("_to_"))
      
      if(nrow(chain1) == nrow(chain2)){
        comb <- mcmc.list(as.mcmc(chain1), as.mcmc(chain2))
        name <- paste0("comb-sim",i,"-b",b,"-c",c)
        
        gelman <- gelman.diag(comb)
        summ <- tibble(true_beta = b, true_clock = c, iteration = i, mpsrf = gelman$mpsrf, 
                       g_clock = gelman$psrf["clock",1], g_beta = gelman$psrf["beta",1],
                       g_01 = gelman$psrf[3,1], g_10 = gelman$psrf[4,1], g_12 = gelman$psrf[5,1], g_21 = gelman$psrf[6,1])
        gelman_sim_thin <- rbind(gelman_sim_thin, summ)
      }
    }
  }
}

View(gelman_sim_thin)



# Compare with true values ----

# calculate means
mean_post <- log_sim_bl1 %>% filter(generation > 10000) %>%
  group_by(true_beta, true_clock, iteration, chain) %>% 
  summarise(mean_clock = mean(clock),
            mean_beta = mean(beta),
            mean_01 = mean(`0_to_1`),
            mean_10 = mean(`1_to_0`),
            mean_12 = mean(`1_to_2`),
            mean_21 = mean(`2_to_1`))

means_chain1 <- filter(mean_post, chain == 1) 
means_chain2 <- filter(mean_post, chain == 2) 

# plot means distribution and mean of means
# clock
ggplot(means_chain1, aes(factor(true_beta), mean_clock)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_clock), col = "red") +
  labs(x = expression("Data simulated under phylogenetic-distance power, " ~ beta), 
       y = expression("Estimated rate of host-repertoire evolution, " ~ mu)) +
  facet_grid(. ~ true_clock) +
  theme_light() +
  theme(strip.text.x = element_text(size = 10, color = "black"),
        strip.background.x = element_rect(fill = "grey90", color = "grey70"))

# beta
ggplot(means_chain1, aes(factor(true_clock), mean_beta)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = true_beta), col = "red") +
  facet_grid(. ~ true_beta) +
  labs(x = expression("Data simulated under host-repertoire evolution rate, " ~ mu), 
       y = expression("Estimated phylogenetic-distance power, " ~ beta)) +
  theme_light() +
  theme(strip.text.x = element_text(size = 10, color = "black"),
        strip.background.x = element_rect(fill = "grey90", color = "grey70")) 

# gain/loss rates
means_rates <- gather(means_chain1,"rate", "value", 7:10) %>% 
  mutate(truth = case_when(rate=="mean_01"~0.05,
                           rate=="mean_10"~0.5,
                           rate=="mean_12"~0.25,
                           rate=="mean_21"~0.2))

#by beta
ggplot(means_rates, aes(factor(true_beta), value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = truth), col = "red") +
  labs(x = expression("Data simulated under phylogenetic-distance power, " ~ beta), 
       y = "Estimated gain/loss rates") +
  facet_grid(. ~ rate) +
  theme_light() +
  theme(strip.text.x = element_text(size = 10, color = "black"),
        strip.background.x = element_rect(fill = "grey90", color = "grey70"))

#by clock
ggplot(means_rates, aes(factor(true_clock), value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = truth), col = "red") +
  labs(x = expression("Data simulated under rate of host-repertoire evolution, " ~ mu), 
       y = "Estimated gain/loss rates") +
  facet_grid(. ~ rate) +
  theme_light() +
  theme(strip.text.x = element_text(size = 10, color = "black"),
        strip.background.x = element_rect(fill = "grey90", color = "grey70"))


# Bayes factor ----

Q <- 0
prior <- rexp(10000, rate = 1)
k_prior <- kdensity(x = prior, kernel='gamma', support=c(0,Inf), bw=0.005)
maximum_a_priori <- k_prior(Q)

BFs <- data.frame()
for(b in betas){
  for(c in clocks){
    for(i in 1:iterations){
      chains <- unique(filter(log_sim_pruned, true_clock == c, true_beta == b, iteration == i)$chain)
      for(ch in chains){
        post <- filter(log_sim_pruned, true_clock == c, true_beta == b, iteration == i, chain == ch)$beta
        k_post <- kdensity(x = post, kernel='gamma', support=c(0,Inf), bw=0.02)
        maximum_a_posteriori = k_post(Q)
        BF <- data.frame(BF = maximum_a_priori/maximum_a_posteriori, true_clock = c, true_beta = b, iteration = i, chain = ch)
        BFs <- rbind(BFs, BF)
      }
    }
  }
}

BFscat <- mutate(BFs, cat = "category")
for(l in 1:nrow(BFscat)){
  bf <- BFscat$BF[l]
  if(is.na(bf)){
    BFscat$cat[l] <- "Decisive"
  } else if(bf <= 1){
    BFscat$cat[l] <- "Favors_M0"
  } else if(bf > 1 & bf <=3){
    BFscat$cat[l] <- "Insubstantial"
  } else if(bf > 3 & bf <= 10){
    BFscat$cat[l] <- "Substantial"
  } else if(bf > 10 & bf <= 30){
    BFscat$cat[l] <- "Strong"
  } else if(bf > 30 & bf <= 100){
    BFscat$cat[l] <- "Very_strong"
  } else if(bf > 100){
    BFscat$cat[l] <- "Decisive"
  }
}

BFscat %>% group_by(true_beta, true_clock, cat) %>% 
  summarize(n = n()) %>% 
  mutate(prop = n/sum(n)) %>% 
  mutate(beta = paste0("beta: ",true_beta)) %>% 
  ggplot(aes(factor(x = true_clock), y = prop*100,
             fill = factor(cat, levels = c("Favors_M0", "Insubstantial", "Substantial", 
                                           "Strong", "Very_strong", "Decisive")))) +
  geom_bar(stat = "identity") +
  scale_fill_grey(labels = function(x) gsub('_', ' ', x)) +
  labs(x = expression("Data simulated under host-repertoire evolution rate, " ~ mu), 
       y = "Percentage of simulations", fill = "Support for MD") +
  facet_grid(. ~ beta, labeller = label_parsed) +
  theme_light() +
  theme(strip.text.x = element_text(size = 10, color = "black"),
        strip.background.x = element_rect(fill = "white", color = "grey70")) 



# Character history ------------------------------------------------------

# function to make a posterior probability matrix
make_anc_node_graph = function(dat, nodes, state) { 
  
  dat <- filter(dat, node_index %in% nodes)
  iterations = sort(unique(dat$iteration))
  n_iter = length(iterations)
  
  # get dimensions
  n_host_tip = length( str_split( dat$start_state[1], "" )[[1]] )
  n_parasite_lineage = length(unique(dat$node_index))
  g = matrix(data = 0, nrow = n_parasite_lineage, ncol = n_host_tip)
  
  for (it in iterations) {
    dat_it = dat[ dat$iteration == it, ]
    ret = list()
    for (i in 1:length(nodes)) {
      dat2 = dat_it[ dat_it$node_index == nodes[i], ]
      if(nrow(dat2)==1){
        ret[[i]] = dat2
      } else{
        ret[[i]] = dat2[ which.min(dat2$transition_time), ]
      }
    }
    
    ret <- rbindlist(ret)
    
    for (r in 1:nrow(ret)) {
      s = as.numeric( str_split (ret$end_state[r], "")[[1]] )
      s_idx = s %in% state
      g[ r, s_idx ] = g[ r, s_idx ] + 1
    }
  }
  
  # convert to probability
  g = g * (1/n_iter)
  
  return(g)
}

path_hist <- "./output/history/"

# list of lists of sampled histories
histories <- tibble() 
burnin <- 1000 # to take the last 1000 rows = 50% burnin
nodes <- 35:67 # select internal nodes

for(b in c(0,1,4)){
  for(c in c(0.1,0.5,1)){
    for(i in 1:50){

      name <- paste0("out.bl1-rep1-sim",i,"-b",b,"-c",c)  # get sampled histories only for chain (rep) 1 for each parameter combination
      
      colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
      history_dat = read.table(paste0(path_hist,name,".history.txt"), sep="\t", header=T, colClasses = colclasses)
      
      its <- distinct(history_dat,iteration) %>% top_n(burnin, iteration) %>% pull()  # select iterations to keep
      history_thin <- filter(history_dat, iteration %in% its) %>%   # thin samples
        mutate(node_index = node_index + 1)   # fix node index - problem with output from RevBayes
      
      pp_s1 <- make_anc_node_graph(history_thin, nodes, c(1)) # get posteriors for state 1
      pp_s2 <- make_anc_node_graph(history_thin, nodes, c(2)) # get posteriors for state 2
      row.names(pp_s1) <- row.names(pp_s2) <- paste0("Index_",nodes)
      colnames(pp_s1) <- colnames(pp_s2) <- hosts
      
      graph1 <- graph_from_incidence_matrix(pp_s1, weighted = TRUE)
      el1 <- get.data.frame(graph1, what = "edges") %>% mutate(p1 = weight) %>% dplyr::select(-weight)
      graph2 <- graph_from_incidence_matrix(pp_s2, weighted = TRUE)
      el2 <- get.data.frame(graph2, what = "edges") %>% mutate(p2 = weight) %>% dplyr::select(-weight)
      
      p_el_full <- full_join(el1,el2, by = c("from","to")) %>% 
        complete(from, to, fill = list(p2 = 0, p1 = 0)) %>%      # cleaner figure but can't see unsampled interactions
        mutate(to = factor(to, levels = hosts_ladder),
               from = factor(from, levels = paste0("Index_",nodes)),
               p3 = p1 + p2, 
               beta = b, clock = c, it = i)
      
      histories <- bind_rows(histories, p_el_full)
    }
  }
}


# Compare truth and posterior ----

mean_infer <- histories %>% 
  group_by(beta, clock, to, from) %>% 
  summarise(m_p0 = 1-mean(p3), m_p1 = mean(p1), m_p2 = mean(p2))

acc_all <- tibble()
for(b in betas){
  for(c in clocks){
    for(i in 1:iterations){
      
      truth <- filter(internal_sims, beta == b, clock == c, it == i)
      infer <- filter(histories, beta == b, clock == c, it == i)
      average <- filter(mean_infer, beta == b, clock == c)
      
      acc <- left_join(truth, infer) %>% 
        mutate(p1 = case_when(is.na(p1)~0,
                              TRUE~p1),
               p2 = case_when(is.na(p2)~0,
                              TRUE~p2),
               p3 = case_when(is.na(p3)~0,
                              TRUE~p3)) %>% 
        left_join(average) %>% 
        mutate(acc = case_when(weight==1 ~ p1,
                               weight==2 ~ p2,
                               weight==0 ~ 1-p3),
               acc_mean = case_when(weight==1 ~ m_p1,
                                    weight==2 ~ m_p2,
                                    weight==0 ~ m_p0),
               beta = b, clock = c, it = i)
      
      acc_all <- bind_rows(acc_all, acc)
    }
  }
}



# Plots

# accuracy
acc_all_plot <- acc_all %>% group_by(beta, clock, it, weight) %>%
  summarise(Observed = mean(acc),
            Expected = mean(acc_mean)) %>%
  gather("group","acc",5:6) %>% 
  ungroup %>% 
  mutate(group = factor(group, levels = c("Observed", "Expected")),
         beta_lab = factor(beta, labels = c(beta ~"= 0", beta ~"= 1", beta ~"= 4")),
         clock_lab = factor(clock, labels = c(mu ~"= 0.1", mu ~"= 0.5", mu ~"= 1.0")))

ggplot(filter(acc_all_plot, weight > 0)) +
  geom_boxplot(aes(weight, acc, col = group)) +
  scale_color_manual(values = c("Expected" = "grey60", "Observed" = "black")) +
  facet_grid(clock_lab ~ beta_lab, labeller = label_parsed) +
  theme_bw() +
  labs(y = "Posterior probability of true ancestral state", x = "True state",
       col = "") +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey40"),
    axis.title.x = element_text(size = 13),
    axis.title.y = element_text(size = 13),
    strip.text.y = element_text(angle = 0))


# Check by comparing truth and inference for one dataset
# truth
filter(internal_sims, beta == 0, clock == 0.1, it == 1) %>% 
  mutate(from = factor(from, levels = hosts_ladder)) %>% 
  ggplot(aes(x = from, y = to, fill = weight)) +
  geom_tile() +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_manual(values = c("grey50", "black")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


# inference
ggplot(filter(histories, beta == 0, clock == 0.1, it == 1)) +
  geom_tile(aes(x = to, y = from, fill = p2)) +
  geom_point(aes(x = to, y = from, col = p1), size = 2, shape = "square") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_color_gradient(low = "white", high = "black") +
  #labs(fill = "Posterior\nprobability") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())

