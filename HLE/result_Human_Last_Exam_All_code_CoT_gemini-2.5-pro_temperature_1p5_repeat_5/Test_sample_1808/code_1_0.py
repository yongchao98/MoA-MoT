import sys

def calculate_fst(p1, p2):
    """Calculates Fst given the allele frequencies in two populations."""
    # Ht: Expected heterozygosity in the total combined population
    p_total = (p1 + p2) / 2.0
    # Avoid division by zero if p_total is 0 or 1
    if p_total == 0.0 or p_total == 1.0:
        return 0.0
    ht = 2.0 * p_total * (1.0 - p_total)
    
    # Hs: Average of expected heterozygosities within each subpopulation
    h1 = 2.0 * p1 * (1.0 - p1)
    h2 = 2.0 * p2 * (1.0 - p2)
    hs = (h1 + h2) / 2.0
    
    # Fst formula
    fst = (ht - hs) / ht
    return fst

def simulate_gene_flow():
    """
    Simulates the effect of gene flow on Fst over time.
    """
    # Initial conditions: Two populations that are completely different.
    # p1 is the frequency of an allele in Population 1.
    # p2 is the frequency of the same allele in Population 2.
    p1 = 0.0
    p2 = 1.0
    
    # Simulation parameters
    migration_rate = 0.05  # 5% of each population migrates each generation
    num_generations = 50
    
    print("This simulation demonstrates that gene flow reduces Fst.")
    print("Initial State: Pop 1 allele freq = 0.0, Pop 2 allele freq = 1.0\n")
    
    for gen in range(num_generations + 1):
        # Calculate Fst for the current generation
        fst = calculate_fst(p1, p2)
        
        if gen % 5 == 0:
            print(f"Generation {gen:2}: Allele Freq (Pop1: {p1:.3f}, Pop2: {p2:.3f}) -> Fst = {fst:.4f}")
        
        # Simulate gene flow for the next generation
        # The new frequency in a population is a weighted average of its
        # current frequency and the frequency of the other population.
        p1_new = p1 * (1 - migration_rate) + p2 * migration_rate
        p2_new = p2 * (1 - migration_rate) + p1 * migration_rate
        
        p1 = p1_new
        p2 = p2_new
        
    print("\nConclusion: As gene flow occurs, Fst drops from 1 (complete differentiation)")
    print("towards 0 (homogenization), showing that high gene flow and high Fst are incompatible.")

if __name__ == '__main__':
    simulate_gene_flow()

<<<A>>>