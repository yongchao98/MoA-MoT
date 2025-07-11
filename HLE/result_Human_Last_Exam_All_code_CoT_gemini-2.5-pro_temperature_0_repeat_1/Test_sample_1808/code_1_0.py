def simulate_gene_flow_on_fst():
    """
    Simulates the effect of gene flow on Fst between two populations.
    """
    # Initial conditions
    # Population 1 is fixed for allele 'A', Population 2 is fixed for allele 'a'
    p1 = 1.0  # Frequency of allele 'A' in population 1
    p2 = 0.0  # Frequency of allele 'A' in population 2
    
    # Migration rate (gene flow)
    m = 0.05  # 5% of each population migrates to the other each generation
    
    # Number of generations to simulate
    generations = [0, 1, 5, 10, 25, 50]
    
    print(f"Simulation of Gene Flow's Effect on Fst")
    print(f"Initial Frequencies: p1 = {p1}, p2 = {p2}")
    print(f"Migration Rate (m) = {m}\n")

    for gen in range(max(generations) + 1):
        # Calculate Fst for the current generation
        # Ht: Expected heterozygosity in the total population
        # Hs: Average expected heterozygosity within subpopulations
        p_bar = (p1 + p2) / 2
        ht = 2 * p_bar * (1 - p_bar)
        hs = (2 * p1 * (1 - p1) + 2 * p2 * (1 - p2)) / 2
        
        # Avoid division by zero if Ht is 0 (populations are fixed for the same allele)
        if ht == 0:
            fst = 0
        else:
            fst = (ht - hs) / ht

        if gen in generations:
            print(f"--- Generation {gen} ---")
            print(f"Allele Frequencies: p1={p1:.3f}, p2={p2:.3f}")
            # The prompt asks to show the numbers in the final equation
            print(f"Fst Calculation: Ht = {ht:.3f}, Hs = {hs:.3f}")
            print(f"Fst = (Ht - Hs) / Ht = ({ht:.3f} - {hs:.3f}) / {ht:.3f} = {fst:.3f}\n")

        # Apply gene flow to update allele frequencies for the next generation
        p1_new = (1 - m) * p1 + m * p2
        p2_new = m * p1 + (1 - m) * p2
        p1, p2 = p1_new, p2_new

simulate_gene_flow_on_fst()