def simulate_gene_flow_and_fst():
    """
    Simulates the effect of gene flow on Fst between two populations.
    """
    # Initial allele 'A' frequencies in population 1 and 2
    p1 = 1.0
    p2 = 0.0

    # Migration rate (gene flow) between populations
    m = 0.05  # 5% of individuals migrate each generation

    # Number of generations to simulate
    generations = 50

    print(f"Initial State:")
    print(f"Allele 'A' Freq in Pop 1 (p1): {p1:.2f}")
    print(f"Allele 'A' Freq in Pop 2 (p2): {p2:.2f}")
    print("-" * 30)
    print("Simulation of Gene Flow:")
    print("Gen | p1   | p2   | Fst")
    print("----|------|------|------")

    for gen in range(generations + 1):
        # Allele 'a' frequencies
        q1 = 1 - p1
        q2 = 1 - p2

        # Calculate Fst = (Ht - Hs) / Ht
        # Hs = average heterozygosity in subpopulations
        hs = (2 * p1 * q1 + 2 * p2 * q2) / 2
        
        # Ht = expected heterozygosity in total population
        p_total = (p1 + p2) / 2
        q_total = (q1 + q2) / 2
        ht = 2 * p_total * q_total

        # Fst calculation
        if ht == 0:
            fst = 0 # Avoid division by zero if populations become identical
        else:
            fst = (ht - hs) / ht

        print(f"{gen:<4}| {p1:<6.3f}| {p2:<6.3f}| {fst:.3f}")

        # Stop if populations are nearly identical
        if fst < 0.001 and gen > 0:
            print("\nPopulations have homogenized. Fst is near zero.")
            break

        # Update allele frequencies due to migration
        p1_new = p1 * (1 - m) + p2 * m
        p2_new = p2 * (1 - m) + p1 * m
        p1 = p1_new
        p2 = p2_new

simulate_gene_flow_and_fst()
<<<A>>>