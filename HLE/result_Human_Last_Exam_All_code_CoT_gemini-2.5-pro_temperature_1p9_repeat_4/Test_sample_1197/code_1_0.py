import numpy as np

def run_simulation(num_genes, is_regional_mutation):
    """
    Simulates Ks and Ka values for a set of genes.
    
    Args:
        num_genes (int): The number of genes to simulate.
        is_regional_mutation (bool): If True, introduces variable mutation rates.

    Returns:
        tuple: A tuple containing lists of ks_values and ka_values.
    """
    # Fixed parameters for the simulation
    synonymous_sites_per_gene = 1000
    nonsynonymous_sites_per_gene = 3000
    
    # Assume purifying selection is constant across genes, reducing Ka relative to the neutral rate
    selection_pressure_coefficient = 0.2 
    
    if is_regional_mutation:
        # SCENARIO B: Regional mutation rates vary from gene to gene (e.g., hotspots/coldspots)
        # This creates an underlying variable that affects both Ks and Ka
        print("Model: Challenging (Regional Mutation Rate)")
        # Each gene gets a different base mutation rate drawn from a distribution
        base_mutation_rates = np.random.uniform(low=0.005, high=0.02, size=num_genes)
    else:
        # STANDARD MODEL: Uniform mutation rate across all genes
        print("Model: Simple Drift (Uniform Mutation Rate)")
        # All genes share the same constant mutation rate
        base_mutation_rates = np.full(num_genes, 0.01)

    # Use the mutation rate to simulate substitutions.
    # The number of substitutions is proportional to the base rate.
    # This is a simplification but captures the core dependency.
    # Ks is a direct reflection of the neutral mutation rate.
    ks_values = base_mutation_rates * (1 + np.random.normal(0, 0.1, num_genes)) # add some noise
    
    # Ka is also affected by the base mutation rate, but suppressed by purifying selection.
    ka_values = (base_mutation_rates * selection_pressure_coefficient) * (1 + np.random.normal(0, 0.1, num_genes)) # add some noise
    
    return ks_values, ka_values

def main():
    """Main function to run simulations and print results."""
    num_simulated_genes = 200

    print("This simulation demonstrates how regional mutation rates can create a correlation")
    print("between synonymous (Ks) and nonsynonymous (Ka) substitution rates,")
    print("which challenges the assumptions of simple genetic drift models.\n")
    
    # --- Run Scenario 1: Uniform Mutation Rate ---
    ks_uniform, ka_uniform = run_simulation(num_simulated_genes, is_regional_mutation=False)
    correlation_uniform = np.corrcoef(ks_uniform, ka_uniform)[0, 1]
    
    print(f"Number of simulated genes = {num_simulated_genes}")
    print("In this model, any correlation is due to random chance, as the underlying mutation rate is constant.")
    print("Final Equation: Correlation(Ks, Ka) = r")
    print(f"Calculated Correlation (r) = {correlation_uniform:.4f}\n")

    print("-" * 50)
    
    # --- Run Scenario 2: Regional (Variable) Mutation Rate ---
    ks_regional, ka_regional = run_simulation(num_simulated_genes, is_regional_mutation=True)
    correlation_regional = np.corrcoef(ks_regional, ka_regional)[0, 1]
    
    print(f"Number of simulated genes = {num_simulated_genes}")
    print("In this model, genes in 'hotspots' have higher Ks and Ka, and genes in 'coldspots' have lower values for both.")
    print("This shared dependency creates a strong positive correlation, illustrating the phenomenon in Answer B.")
    print("Final Equation: Correlation(Ks, Ka) = r")
    print(f"Calculated Correlation (r) = {correlation_regional:.4f}\n")

    print("The high correlation in the second scenario makes it difficult to parse the effects of")
    print("selection from the effects of mutational bias, thus challenging the models.")

if __name__ == "__main__":
    main()