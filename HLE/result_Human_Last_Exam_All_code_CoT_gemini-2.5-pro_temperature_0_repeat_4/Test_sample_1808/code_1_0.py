import numpy as np

def calculate_fst(p1, p2):
    """Calculates Fst for two populations with given allele frequencies."""
    p_bar = (p1 + p2) / 2.0
    # Avoid division by zero if p_bar is 0 or 1 (no variation in total pop)
    if p_bar == 0 or p_bar == 1:
        return 0.0
    var_p = ((p1 - p_bar)**2 + (p2 - p_bar)**2) / 2.0
    fst = var_p / (p_bar * (1 - p_bar))
    return fst

def simulate_gene_flow():
    """
    Simulates the effect of gene flow on Fst over generations.
    """
    # Initial state:
    # Population 1 is fixed for allele 'A'
    p1 = 1.0
    # Population 2 is fixed for allele 'a' (frequency of 'A' is 0)
    p2 = 0.0
    
    # Migration rate (m): 5% of each population migrates to the other each generation
    m = 0.05
    
    num_generations = 20

    print("This simulation demonstrates the effect of gene flow on Fst.")
    print(f"Initial state: Pop 1 freq(A) = {p1}, Pop 2 freq(A) = {p2}")
    print(f"Migration rate (m) = {m}\n")

    for gen in range(num_generations + 1):
        # Calculate Fst for the current generation
        fst = calculate_fst(p1, p2)
        print(f"Generation {gen:2d}: Fst = {fst:.4f}")
        
        # Stop if we are at the last generation
        if gen == num_generations:
            break

        # Apply gene flow to calculate allele frequencies for the next generation
        # The new frequency in pop 1 is a mix of the non-migrants from pop 1 and migrants from pop 2
        p1_next = (1 - m) * p1 + m * p2
        # The new frequency in pop 2 is a mix of the non-migrants from pop 2 and migrants from pop 1
        p2_next = (1 - m) * p2 + m * p1
        
        p1, p2 = p1_next, p2_next
        
    print("\nConclusion: As gene flow occurs over generations, Fst decreases.")
    print("Therefore, high gene flow is incompatible with a high Fst.")

if __name__ == "__main__":
    simulate_gene_flow()
