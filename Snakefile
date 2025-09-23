SEEDS    = [9,20,31,42,53,64,75,86,97,108,119,130,141,152,163]   # 15 seeds
CORRUPTS = [10, 50, 100]
RATES    = [0.3, 0.4, 0.5, 0.6, 0.7]
BETAS  = [0, 1, 4]
CLOCKS = [0.1, 0.5, 1]



HOST_TREE  = "angio_25tips_origin.phy"
VIRUS_TREE = "Nymphalini.phy"

# --------------------
# Final targets
# --------------------
rule all:
    input:
        "precision_by_corrupt.png",
        "recall_by_corrupt.png",
        "f1_by_corrupt.png",
        "precision_by_corrupt_elbow.png",
        "recall_by_corrupt_elbow.png",
        "f1_by_corrupt_elbow.png",
        expand(
            "results/braga/{seed}/b{beta}_c{clock}/metrics.json",
            seed=SEEDS, beta=BETAS, clock=CLOCKS
        ),
        "braga_precision.png",
        "braga_recall.png",
        "braga_f1.png"

# --------------------
# Simulation rule
# --------------------
rule simulate:
    output:
        original  = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/original_matrix.png",
        flipped   = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/flipped_matrix.png",
        metrics   = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/metrics.json",
        elbow_metrics = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/elbow_metrics.json",
        log_out   = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/simulate.out",
        log_err   = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/simulate.err"
    params:
        outdir = "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}"
    shell:
        """
        mkdir -p {params.outdir}
        {{
            python simulation.py \
                --seed {wildcards.seed} \
                --host_tree {HOST_TREE} \
                --virus_tree {VIRUS_TREE} \
                --corrupt {wildcards.corrupt} \
                --r01_p {wildcards.r01} \
                --r10_p {wildcards.r10} \
                --r01_h {wildcards.r01} \
                --r10_h {wildcards.r10} \
                --outdir {params.outdir} \
                --original_metrics {output.metrics} \
                --elbow_metrics {output.elbow_metrics} \
                > {output.log_out} 2> {output.log_err}
        }}
        """


rule braga_simulator:
    output:
        nexus    = "results/braga/{seed}/sim{seed}-b{beta}-c{clock}.nex",
        settings = "results/braga/{seed}/sim{seed}-b{beta}-c{clock}.settings.txt"
    params:
        outdir = "results/braga/{seed}"
    shell:
        """
        mkdir -p {params.outdir}
        rb Simulate.Rev \
            --args {wildcards.seed} {wildcards.beta} {wildcards.clock} {params.outdir} \
            > {params.outdir}/braga_{wildcards.seed}_b{wildcards.beta}_c{wildcards.clock}.log 2>&1
        """



rule braga_pipeline:
    input:
        nexus    = "results/braga/{seed}/sim{seed}-b{beta}-c{clock}.nex",
        settings = "results/braga/{seed}/sim{seed}-b{beta}-c{clock}.settings.txt"
    output:
        original = "results/braga/{seed}/b{beta}_c{clock}/original_matrix.png",
        flipped  = "results/braga/{seed}/b{beta}_c{clock}/flipped_matrix.png",
        metrics  = "results/braga/{seed}/b{beta}_c{clock}/metrics.json",
        log_out  = "results/braga/{seed}/b{beta}_c{clock}/simulate.out",
        log_err  = "results/braga/{seed}/b{beta}_c{clock}/simulate.err"
    params:
        outdir = "results/braga/{seed}/b{beta}_c{clock}"
    shell:
        """
        mkdir -p {params.outdir}
        python simulation.py \
            --seed {wildcards.seed} \
            --host_tree {HOST_TREE} \
            --virus_tree {VIRUS_TREE} \
            --nexus_file {input.nexus} \
            --braga {input.settings} \
            --outdir {params.outdir} \
            --original_metrics {output.metrics} \
            > {output.log_out} 2> {output.log_err}
        """

# --------------------
# Plotting rule
# --------------------
rule plot_my_simulation:
    input:
        expand(
            "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/metrics.json",
            seed=SEEDS, corrupt=CORRUPTS, r01=RATES, r10=RATES
        ),
        expand(
            "results/{seed}/corrupt{corrupt}/r01{r01}_r10{r10}/elbow_metrics.json",
            seed=SEEDS, corrupt=CORRUPTS, r01=RATES, r10=RATES
        )
    output:
        "precision_by_corrupt.png",
        "recall_by_corrupt.png",
        "f1_by_corrupt.png",
        "precision_by_corrupt_elbow.png",
        "recall_by_corrupt_elbow.png",
        "f1_by_corrupt_elbow.png"
    shell:
        """
        python plotting.py
        """
rule plot_braga_simulation:
    input:
        expand(
            "results/braga/{seed}/b{beta}_c{clock}/metrics.json",
            seed=SEEDS, beta=BETAS, clock=CLOCKS
        )
    output:
        "braga_precision.png",
        "braga_recall.png",
        "braga_f1.png",
    shell:
        """
        python braga_plotting.py
        """
