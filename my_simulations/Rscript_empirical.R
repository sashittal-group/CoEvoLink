# Script for analysing RevBayes output for empirical study performed in Braga et al. 2020
# Bayesian inference of ancestral host-parasite interactions under a phylogenetic model of host repertoire evolution

library(ape)
library(tidyverse)
library(igraph)
library(MCMCpack)
library(coda)
library(kdensity)
library(phylotate)
library(data.table)

# Read trees ----
tree <- read.tree("Nymphalini.phy")
host_tree <- read.tree("angio_25tips_bl1.phy")
nhosts <- Ntip(host_tree)
hosts <- host_tree$tip.label

# Convergence ----
# 3 state dataset
# Host tree with all branch lenghts set to 1
chain1 <- read.table("out.2.bl1.nymphalini.3s.log", header = TRUE)[,c(1,5,7:11)]
chain2 <- read.table("out.9.bl1.nymphalini.3s.log", header = TRUE)[,c(1,5,7:11)]
gelman.diag(mcmc.list(as.mcmc(filter(chain1, Iteration >= 100000 & Iteration <= 500000)), 
                      as.mcmc(filter(chain2, Iteration >= 100000 & Iteration <= 500000))))

gelman.plot(mcmc.list(as.mcmc(filter(chain1, Iteration >= 100000 & Iteration <= 500000)[,-1]),
                      as.mcmc(filter(chain2, Iteration >= 100000 & Iteration <= 500000)[,-1])))

# Host tree with original branch lenghts (time calibrated)
chain3 <- read.table("out.2.real.nymphalini.3s.log", header = TRUE)[,c(1,5,7:11)]
chain4 <- read.table("out.9.real.nymphalini.3s.log", header = TRUE)[,c(1,5,7:11)]
gelman.diag(mcmc.list(as.mcmc(filter(chain3, Iteration >= 100000 & Iteration <= 500000)), 
                      as.mcmc(filter(chain4, Iteration >= 100000 & Iteration <= 500000))))

gelman.plot(mcmc.list(as.mcmc(filter(chain3, Iteration >= 100000 & Iteration <= 500000)[,-1]),
                      as.mcmc(filter(chain4, Iteration >= 100000 & Iteration <= 500000)[,-1])))


# 2 state dataset
# Host tree with all branch lenghts set to 1
chain5 <- read.table("out.4.bl1.nymphalini.2s.log", header = TRUE)[,c(1,5,7:11)]
chain6 <- read.table("out.9.bl1.nymphalini.2s.log", header = TRUE)[,c(1,5,7:11)]
gelman.diag(mcmc.list(as.mcmc(filter(chain5, Iteration <= 100000)), 
                      as.mcmc(filter(chain6, Iteration <= 100000))))

gelman.plot(mcmc.list(as.mcmc(filter(chain5, Iteration <= 100000)[,-1]),
                      as.mcmc(filter(chain6, Iteration <= 100000)[,-1])))

# Host tree with original branch lenghts (time calibrated)
chain7 <- read.table("out.4.real.nymphalini.2s.log", header = TRUE)[,c(1,5,7:11)]
chain8 <- read.table("out.9.real.nymphalini.2s.log", header = TRUE)[,c(1,5,7:11)]
gelman.diag(mcmc.list(as.mcmc(filter(chain7, Iteration <= 96500)), 
                      as.mcmc(filter(chain8, Iteration <= 96500))))

gelman.plot(mcmc.list(as.mcmc(filter(chain7, Iteration <= 96500)[,-1]),
                      as.mcmc(filter(chain8, Iteration <= 96500)[,-1])))


post3sb <- filter(chain1, Iteration >= 100000 & Iteration <= 500000)
post3sr <- filter(chain3, Iteration >= 100000 & Iteration <= 500000)
post2sb <- filter(chain5, Iteration >= 100000 & Iteration <= 500000)
post2sr <- filter(chain7, Iteration >= 100000 & Iteration <= 500000)

colnames(post3sb) <- colnames(post3sr) <- c("generation","clock","beta", "lambda[01]", "lambda[10]", "lambda[12]", "lambda[21]")
colnames(post2sb) <- colnames(post2sr) <- c("generation","clock","beta", "lambda[01]", "lambda[10]", "lambda[12]", "lambda[21]")

posterior <- bind_rows(post3sb,post3sr,post2sb,post2sr) %>% 
  mutate(states = c(rep(3,16002),rep(2,16002)),
         tree = c(rep("bl1", 8001),rep("real", 8001),rep("bl1", 8001),rep("real", 8001)))


# Parameter estimates ----

# mean estimates
means <- group_by(posterior, states, tree) %>% 
  summarise_all(mean)

# Density

prior_beta <- data.frame(sim = NA, beta = rexp(10000, rate = 1))
prior_clock <- data.frame(sim = NA, clock = rexp(10000, rate = 10))

# get density function for the other parameters

for(m in c("real", "bl1")){
  for(s in 2:3){
    for(p in colnames(posterior)[c(2:7)]){
      post <- filter(posterior, states == s, tree == m) %>% dplyr::select(p)
      if(p %in% c("beta", "clock")) {
        kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.01)
      } else{
        kdens <- kdensity(x = post[,1], kernel='gamma', support=c(0,Inf), bw = 0.001)
      }
      assign(paste0("kd_",p,"_",m,"_",s,"s"),kdens)
    }
  }
}

# create tibbles and plot
z = seq(0,5,0.001)
y = seq(0,1,0.0002)

# beta
dens_beta <- tibble(x = z, prior = dexp(x = z, rate=1), real_3s = kd_beta_real_3s(x), real_2s = kd_beta_real_2s(x), 
                    bl1_3s = kd_beta_bl1_3s(x), bl1_2s = kd_beta_bl1_2s(x))

ggplot(dens_beta) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  labs(x = expression("Estimated phylogenetic-distance power, " ~ beta)) +
  theme_light()

# clock
dens_clock <- tibble(x = z, prior = dexp(x = z, rate=10), real_3s = kd_clock_real_3s(x), real_2s = kd_clock_real_2s(x), 
                     bl1_3s = kd_clock_bl1_3s(x), bl1_2s = kd_clock_bl1_2s(x))
ggplot(dens_clock) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  scale_x_continuous(limits = c(0,3)) +
  labs(x = expression("Estimated rate of host-repertoire evolution, " ~ mu)) +
  theme_light()

# rates
prior_rates <- dbeta(y, 1, 3)
dens_01 <- tibble(x = y, prior = prior_rates, real_3s = `kd_lambda[01]_real_3s`(x), real_2s = `kd_lambda[01]_real_2s`(x), 
                  bl1_3s = `kd_lambda[01]_bl1_3s`(x), bl1_2s = `kd_lambda[01]_bl1_2s`(x))
ggplot(dens_01) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  scale_x_continuous(limits = c(0,0.05)) +
  theme_light()

dens_10 <- tibble(x = y, prior = prior_rates, real_3s = `kd_lambda[10]_real_3s`(x), real_2s = `kd_lambda[10]_real_2s`(x), 
                  bl1_3s = `kd_lambda[10]_bl1_3s`(x), bl1_2s = `kd_lambda[10]_bl1_2s`(x))
ggplot(dens_10) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_light()

dens_12 <- tibble(x = y, prior = prior_rates, real_3s = `kd_lambda[12]_real_3s`(x), real_2s = `kd_lambda[12]_real_2s`(x), 
                  bl1_3s = `kd_lambda[12]_bl1_3s`(x), bl1_2s = `kd_lambda[12]_bl1_2s`(x))
ggplot(dens_12) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  scale_x_continuous(limits = c(0,1)) +
  theme_light()

dens_21 <- tibble(x = y, prior = prior_rates, real_3s = `kd_lambda[21]_real_3s`(x), real_2s = `kd_lambda[21]_real_2s`(x), 
                  bl1_3s = `kd_lambda[21]_bl1_3s`(x), bl1_2s = `kd_lambda[21]_bl1_2s`(x))
ggplot(dens_21) +
  geom_line(aes(x, prior), col = "grey50", alpha = 0.8) +
  geom_line(aes(x, bl1_3s), col = "blue", alpha = 0.8) +
  geom_line(aes(x, bl1_2s), col = "orange", alpha = 0.8) +
  scale_x_continuous(limits = c(0,0.4)) +
  theme_light()


# Bayes factor ----

d_prior <- dexp(x=0, rate=1)

max_3s_bl1 = kd_beta_bl1_3s(0)
max_3s_real = kd_beta_real_3s(0)

BF3s_bl1 <- d_prior/max_3s_bl1
BF3s_real <- d_prior/max_3s_real

max_2s_bl1 = kd_beta_bl1_2s(0)
max_2s_real = kd_beta_real_2s(0)

BF2s_bl1 <- d_prior/max_2s_bl1
BF2s_real <- d_prior/max_2s_real


# Character histories ----

# function to make posterior probability matrix
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

name <- "out.2.bl1.nymphalini.3s.history.txt"
name <- "out.2.real.nymphalini.3s.history.txt"
name <- "out.4.bl1.nymphalini.2s.history.txt"
name <- "out.4.real.nymphalini.2s.history.txt"

colclasses <- c(rep("numeric",7),"character","character","numeric","character",rep("numeric",3))
history_dat = read.table(name, sep="\t", header=T, colClasses = colclasses)

history_dat <- filter(history_dat, iteration >= 100000 & iteration <= 500000) %>% 
  mutate(node_index = node_index + 1)


# Calculate effective rate of evolution ----
tree_length <- sum(tree$edge.length)+20
n_events <- group_by(history_dat,iteration) %>% 
  summarise(n = n()) %>% 
  summarise(mean = mean(n)) %>% 
  pull(mean)

rate <- n_events/tree_length
rate
#---

nodes <- 35:67 # select all internal nodes

# Calculate posterior probability for ancestral states - separately for states 1 and 2
pp_s1 <- make_anc_node_graph(history_dat, nodes, c(1))
row.names(pp_s1) <- paste0("Index_",nodes)
colnames(pp_s1) <- hosts
assign(paste0("pp_s1_",name), pp_s1)

pp_s2 <- make_anc_node_graph(history_dat, nodes, c(2))
row.names(pp_s2) <- paste0("Index_",nodes)
colnames(pp_s2) <- hosts
assign(paste0("pp_s2_",name), pp_s2)

graph1 <- graph_from_incidence_matrix(pp_s1, weighted = TRUE)
el1 <- get.data.frame(graph1, what = "edges") %>% mutate(p1 = weight) %>% dplyr::select(-weight)
graph2 <- graph_from_incidence_matrix(pp_s2, weighted = TRUE)
el2 <- get.data.frame(graph2, what = "edges") %>% mutate(p2 = weight) %>% dplyr::select(-weight)

# Combine states 1 and 2 and calculate p3, the probability of being on fundamental host repertoire
p_el_full <- full_join(el1,el2, by = c("from","to")) %>% 
  complete(nesting(from,to), fill = list(p2 = 0)) %>%      # cleaner figure but can't see unsampled interactions
  mutate(to = factor(to, levels = hosts),
         from = factor(from, levels = paste0("Index_",nodes)),
         p3 = p1+p2)

ggplot() +
  geom_tile(aes(x = to, y = from, fill = p2), data = filter(p_el_full, p2 >= 0)) + # adjust probability threshold
  geom_point(aes(x = to, y = from, col = p3), data = filter(p_el_full, p3 >= 0), size = 2, shape = "square") +
  scale_x_discrete(drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  scale_fill_gradient(low = "white", high = "black") +
  scale_color_gradient(low = "white", high = "black") +
  labs(fill = "Posterior\nprobability") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 270, hjust = 0, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank())


