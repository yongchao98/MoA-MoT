def calculate_fst(N, m):
    """
    Calculates the expected Fst value under the island model of migration.
    
    Args:
        N (int): Effective population size.
        m (float): Migration rate (proportion of a population made up of migrants each generation).
        
    Returns:
        float: The calculated Fst value.
    """
    # The term 4Nm represents the effective number of migrants per generation.
    # If m is 0 (no gene flow), Fst is 1 (complete differentiation).
    # As m increases, Fst approaches 0.
    if 4 * N * m < 0:
        return None # Invalid input
    return 1 / (1 + 4 * N * m)

def demonstrate_gene_flow_effect_on_fst():
    """
    Prints a table showing how Fst changes with different gene flow rates.
    """
    # Assume a constant effective population size
    N = 1000 
    
    # Define a range of migration rates (m) to represent different levels of gene flow
    # m = 0.0001 -> Very low gene flow
    # m = 0.01   -> Moderate gene flow
    # m = 0.1    -> High gene flow
    migration_rates = [0, 0.0001, 0.001, 0.01, 0.1]
    
    print(f"Demonstrating the effect of gene flow on Fst (assuming N={N}):")
    print("-" * 45)
    print(f"{'Migration Rate (m)':<20} | {'Fst (Differentiation)':<25}")
    print("-" * 45)
    
    for m in migration_rates:
        fst = calculate_fst(N, m)
        description = ""
        if m == 0:
            description = "No gene flow"
        elif m <= 0.001:
            description = "Low gene flow"
        elif m <= 0.01:
            description = "Moderate gene flow"
        else:
            description = "High gene flow"
            
        print(f"{m:<20.4f} | {fst:<25.4f}")

    print("-" * 45)
    print("\nAs the migration rate (m) increases, Fst decreases.")
    print("This shows that high gene flow leads to low population differentiation (Fst).")
    print("Therefore, high Fst cannot occur when there is significant gene flow.")

if __name__ == '__main__':
    demonstrate_gene_flow_effect_on_fst()