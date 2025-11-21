# CoEvoLink

A computational framework for analyzing host-parasite/virus coevolutionary relationships using phylogenetic networks and graph-theoretic optimization.

## Overview

CoEvoLink implements novel algorithms for inferring host-virus interaction patterns by combining phylogenetic information from both host and parasite/virus lineages. The framework uses network cut algorithms to identify and recover missing or uncertain interactions in host-parasite association matrices.

### Key Features

- **Phylogenetic Network Analysis**: Integrates phylogenetic trees from both hosts and parasites/viruses
- **Graph Cut Optimization**: Uses minimum cut algorithms to infer missing interactions
- **Parsimony-Based Inference**: Employs Sankoff parsimony for phylogenetic character reconstruction
- **Simulation Framework**: Comprehensive simulation tools for validating methods
- **Multiple Real-World Datasets**: Includes analysis on virion, CLOVER, and other databases

## Installation

### Prerequisites

- Python 3.8+ (tested with Python 3.8-3.12)
- Snakemake (for pipeline execution)
- RevBayes (for Bayesian phylogenetic simulations)

**Note**: This software has been primarily developed and tested with Python 3.8-3.12. While newer versions should work, they have not been extensively tested.

### Python Dependencies

The project requires the following Python packages:

```bash
pip install numpy scipy matplotlib biopython networkx pandas kneed
```

Key dependencies:
- `biopython`: Phylogenetic tree manipulation
- `networkx`: Graph algorithms and minimum cut
- `numpy`, `scipy`: Numerical computations
- `matplotlib`: Visualization
- `kneed`: Elbow detection for parameter optimization

**Note**: For reproducible installations, consider creating a virtual environment. A `requirements.txt` file may be added in future releases for version-pinned dependencies.

## Repository Structure

```
.
├── simulation_utils.py      # Core simulation and analysis utilities
├── Snakefile               # Main Snakemake workflow
├── my_simulations/         # Simulation experiments and results
├── my_simulations_bigger_tree/  # Large-scale simulation tests
├── virion/                 # Virion database analysis
├── clover/                 # CLOVER database analysis
├── aPhyloGeo/             # aPhyloGeo analysis
├── Fresnel/               # Fresnel ensemble methods
├── HostPredictionReview/  # Host prediction comparison studies
├── virion_clover_overlap/ # Overlap analysis between databases
└── Supplementary.pdf      # Supplementary materials
```

## Usage

### Running Simulations

#### Basic Simulation

```bash
cd my_simulations
python simulation.py \
    --seed 42 \
    --host_tree angio_25tips_origin.phy \
    --virus_tree Nymphalini.phy \
    --corrupt 100 \
    --r01_p 0.5 \
    --r10_p 0.5 \
    --r01_h 0.5 \
    --r10_h 0.5 \
    --outdir results/
```

Parameters:
- `--seed`: Random seed for reproducibility
- `--host_tree`: Newick format host phylogenetic tree
- `--virus_tree`: Newick format virus/parasite phylogenetic tree
- `--corrupt`: Number of interactions to randomly hide
- `--r01_p`, `--r10_p`: Transition rates (0→1, 1→0) for parasites
- `--r01_h`, `--r10_h`: Transition rates for hosts
- `--outdir`: Output directory for results

#### Using Snakemake Pipeline

Run the complete analysis pipeline:

```bash
snakemake --cores 4
```

This will:
1. Generate simulated host-parasite networks
2. Corrupt interaction matrices
3. Recover interactions using the network cut algorithm
4. Evaluate performance metrics (precision, recall, F1)
5. Generate visualization plots

### Working with Real Data

#### Virion Database Analysis

```bash
cd virion/
snakemake --cores 4
```

Analyzes virus-host associations from the Virion database.

#### CLOVER Database Analysis

```bash
cd clover
python simulation.py \
    --seed 42 \
    --host_tree <host_tree_file.nwk> \
    --virus_tree <virus_tree_file.nwk> \
    --association_csv CLOVER_0.1_MammalViruses_AssociationsFlatFile.csv \
    --outdir clover_bat_cov_output
```

Processes mammal-virus associations from the CLOVER database. The script requires:
- Host phylogenetic tree in Newick format
- Virus phylogenetic tree in Newick format
- CSV file with Host-Virus associations


## Methodology

### Network Cut Algorithm

The core algorithm constructs a flow network where:

1. **Nodes**: Internal nodes from host and parasite phylogenies + leaf interaction cells
2. **Edges**: Weighted by phylogenetic similarity (transition probabilities)
3. **Source/Sink**: Connected to cells based on observed interaction state
4. **Objective**: Minimize cut value (balance phylogenetic concordance and data fit)

The algorithm finds an optimal balance between:
- **Phylogenetic concordance**: Minimizing changes along branches (parameter λ)
- **Data fidelity**: Maintaining observed interactions (parameter 1-λ)

### Elbow Detection

The framework uses automatic elbow detection (via `kneed` package) to select the optimal λ parameter by analyzing the parsimony-vs-flips curve.

## Output Files

Typical output includes:

- `original_matrix.png`: Visualization of true interaction matrix
- `flipped_matrix.png`: Recovered interaction matrix with highlights
- `metrics.json`: Performance metrics (precision, recall, F1)
- `*.out`, `*.err`: Execution logs
- `*_by_corrupt.png`: Performance plots across corruption levels

## Example Datasets

The repository includes several example datasets:

1. **Angiosperm-Nymphalini**: Plant-butterfly associations
   - Host: `angio_25tips_origin.phy` (25 angiosperm taxa)
   - Virus: `Nymphalini.phy` (butterfly taxa)

2. **Virion**: Comprehensive virus-host database
   - Multiple virus families
   - Bat and human hosts

3. **CLOVER**: Mammal-virus associations
   - Large-scale dataset
   - Multiple host and virus taxonomic groups

## Performance Evaluation

The framework evaluates inference quality using:

- **Precision**: Fraction of predicted interactions that are correct
- **Recall**: Fraction of true interactions that are recovered
- **F1 Score**: Harmonic mean of precision and recall
- **Parsimony Score**: Total evolutionary changes across both phylogenies

## Citation

If you use CoEvoLink in your research, please cite:

```
[Citation information will be added upon publication. Preprint/manuscript in preparation.]
```

For now, you can reference this repository:
```
CoEvoLink: Phylogenetic Network Analysis for Host-Parasite Coevolution
https://github.com/sashittal-group/CoEvoLink
```

## Comparison with Other Methods

The repository includes comparisons with:
- **WISH**: Host prediction method
- **PHiST**: Phylogenetic similarity method
- **PBLks**: Pattern-based approach
- **Virion Ensemble**: Machine learning ensemble
- **ParaFit/PACo**: Cophylogenetic methods (via aPhyloGeo)

Results are available in the `HostPredictionReview/` directory.

## Advanced Features

### Braga Simulator Integration

For more realistic simulations using RevBayes:

```bash
rb Simulate.Rev --args <seed> <beta> <clock> <outdir>
```

### Custom Weight Matrices

You can provide custom weight matrices for edges in the network cut algorithm by modifying the `build_weight_matrix()` function in `simulation_utils.py`.

### Parameter Sweeps

Run parameter sweeps across multiple corruption levels and evolutionary rates:

```bash
snakemake --config corrupt=10,50,100 rates=0.3,0.5,0.7
```

## Troubleshooting

### Common Issues

1. **Missing Dependencies**: Install all required packages via pip
2. **Tree Format Errors**: Ensure phylogenetic trees are in valid Newick format
3. **Memory Issues**: For large trees (>1000 tips), consider using the `my_simulations_bigger_tree/` optimized version
4. **Snakemake Errors**: Check that all input files exist in the expected locations

## Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch
3. Make your changes with clear commit messages
4. Submit a pull request

## License

This project is currently under development. License information will be specified by the Sashittal Lab. 

For usage inquiries, please contact the repository maintainers or open an issue on GitHub.

## Contact

For questions, issues, or collaboration inquiries, please open an issue on GitHub or contact the Sashittal Lab.

## Acknowledgments

This work uses data from:
- **Virion**: Viral host prediction database
- **CLOVER**: Mammal virus associations database
- **aPhyloGeo**: Geographic and phylogenetic analysis tools

Special thanks to all contributors and collaborators in the development of these methods.
