import numpy as np

def simulate_and_explain_correlation():
    """
    This script illustrates why the correlation between synonymous (dS) and
    nonsynonymous (dN) substitution rates challenges simple genetic drift models.

    The core idea is that mutation rates can vary across the genome. Regions
    with high mutation rates ("hotspots") will accumulate both synonymous and
    nonsynonymous mutations faster than regions with low rates ("coldspots").
    This creates a correlation between dS and dN that is not due to random drift
    alone, but to this underlying regional variation.
    """
    print("--- Simulating Gene Evolution in Different Genomic Regions ---")

    # --- Parameters ---
    # Assume nonsynonymous sites are under purifying selection, so their
    # substitution rate is a fraction of the base mutation rate.
    purifying_selection_factor = 0.2  # dN is 20% of the base mutation rate

    # Define base mutation rates for different genomic regions
    coldspot_base_rate = 0.01
    hotspot_base_rate = 0.05

    # --- Calculations ---
    # In this simple model, we assume dS is a direct proxy for the neutral mutation rate.
    # dN is the neutral mutation rate modulated by selection.

    # For a gene in a "coldspot"
    dS_coldspot = coldspot_base_rate
    dN_coldspot = coldspot_base_rate * purifying_selection_factor

    # For a gene in a "hotspot"
    dS_hotspot = hotspot_base_rate
    dN_hotspot = hotspot_base_rate * purifying_selection_factor

    print("\n1. Individual Gene Calculations:")
    print(f"Gene in Coldspot (low mutation rate = {coldspot_base_rate}):")
    print(f"  dS = {dS_coldspot:.2f}")
    print(f"  dN = (base rate) * (selection factor) => {coldspot_base_rate:.2f} * {purifying_selection_factor:.2f} = {dN_coldspot:.3f}")

    print(f"\nGene in Hotspot (high mutation rate = {hotspot_base_rate}):")
    print(f"  dS = {dS_hotspot:.2f}")
    print(f"  dN = (base rate) * (selection factor) => {hotspot_base_rate:.2f} * {purifying_selection_factor:.2f} = {dN_hotspot:.3f}")

    # --- Demonstrating the Correlation ---
    # To statistically show the correlation, we simulate a larger set of genes.
    print("\n2. Demonstrating Correlation Across Many Genes:")
    num_genes = 50
    # Generate random base mutation rates for each gene's region
    np.random.seed(0)
    base_rates = np.random.uniform(0.005, 0.06, num_genes)
    
    # Calculate dS and dN for all genes
    dS_values = base_rates
    dN_values = base_rates * purifying_selection_factor
    
    # Calculate the correlation coefficient
    correlation_matrix = np.corrcoef(dS_values, dN_values)
    correlation = correlation_matrix[0, 1]

    print(f"Simulated {num_genes} genes with varying regional mutation rates.")
    print(f"The calculated correlation coefficient between dS and dN is: {correlation:.4f}")
    print("\nConclusion: A strong positive correlation is observed. This shows that regional mutation variation, a non-random factor, can systematically link dS and dN, challenging predictive models based on simple, uniform genetic drift.")

# Run the simulation
simulate_and_explain_correlation()