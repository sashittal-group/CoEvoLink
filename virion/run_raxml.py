#!/usr/bin/env python3
"""
Script to run RAxML-NG on FASTA sequence files
"""
import os
import sys
import subprocess

def run_raxml_ng(input_fasta, output_prefix, model="GTR+G", bootstrap=1000, threads="auto"):
    """
    Run RAxML-NG on a FASTA alignment file
    
    Parameters:
    input_fasta: Path to input aligned FASTA file
    output_prefix: Prefix for output files
    model: Evolutionary model (GTR+G for DNA, PROTGTR+G for proteins)
    bootstrap: Number of bootstrap replicates
    threads: Number of threads to use
    """
    
    # Construct RAxML-NG command
    cmd = [
        "raxml-ng",
        "--all",
        "--msa", input_fasta,
        "--model", model,
        "--prefix", output_prefix,
        "--bs-trees", str(bootstrap),
        "--threads", str(threads),
        "--seed", "12345"
    ]
    
    print(f"Running RAxML-NG command:")
    print(" ".join(cmd))
    print()
    
    try:
        # Run the command
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print("RAxML-NG completed successfully!")
        print("Output files created:")
        print(f"  - {output_prefix}.raxml.bestTree (best ML tree)")
        print(f"  - {output_prefix}.raxml.support (tree with bootstrap support)")
        print(f"  - {output_prefix}.raxml.bestModel (best-fit model)")
        print(f"  - {output_prefix}.raxml.log (log file)")
        return True
        
    except subprocess.CalledProcessError as e:
        print(f"Error running RAxML-NG: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return False
    except FileNotFoundError:
        print("Error: raxml-ng not found. Please install RAxML-NG first.")
        print("Install with: conda install -c bioconda raxml-ng")
        return False

def main():
    # Define input files and output prefixes
    alignments_dir = "virion_bat_cov_output/alignments"
    
    analyses = [
        {
            "input": f"{alignments_dir}/hosts_aligned.fasta",
            "prefix": "virion_bat_cov_output/trees/hosts_raxml",
            "model": "GTR+G"  # DNA model for host sequences
        },
        {
            "input": f"{alignments_dir}/viruses_aligned.fasta", 
            "prefix": "virion_bat_cov_output/trees/viruses_raxml",
            "model": "PROTGTR+G"  # Protein model for virus sequences
        }
    ]
    
    # Create output directory if it doesn't exist
    os.makedirs("virion_bat_cov_output/trees", exist_ok=True)
    
    # Run RAxML-NG for each alignment
    for analysis in analyses:
        if os.path.exists(analysis["input"]):
            print(f"\n{'='*60}")
            print(f"Processing: {analysis['input']}")
            print(f"Output prefix: {analysis['prefix']}")
            print(f"Model: {analysis['model']}")
            print(f"{'='*60}")
            
            success = run_raxml_ng(
                input_fasta=analysis["input"],
                output_prefix=analysis["prefix"],
                model=analysis["model"],
                bootstrap=1000,
                threads="auto"
            )
            
            if not success:
                print(f"Failed to process {analysis['input']}")
        else:
            print(f"Warning: {analysis['input']} not found, skipping...")

if __name__ == "__main__":
    main()