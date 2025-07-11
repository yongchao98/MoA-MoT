import random

def simulate_population_genetics():
    """
    Simulates and compares genetic differentiation for autosomal vs. Y-linked markers
    between males and females in a hypothetical population.
    """
    n_males = 500
    n_females = 500

    # --- Case 1: Autosomal Marker ---
    # This marker is on a non-sex chromosome and is inherited by both sexes.
    # Allele frequencies in the overall gene pool are assumed to be equal.
    # We expect allele frequencies in male and female subgroups to be very similar.
    p_autosomal = 0.4  # Frequency of allele 'T'
    
    male_T_alleles = 0
    female_T_alleles = 0

    # Simulate genotypes for males (2 alleles each)
    for _ in range(n_males):
        # Draw two alleles from the gene pool
        if random.random() < p_autosomal:
            male_T_alleles += 1
        if random.random() < p_autosomal:
            male_T_alleles += 1

    # Simulate genotypes for females (2 alleles each)
    for _ in range(n_females):
        if random.random() < p_autosomal:
            female_T_alleles += 1
        if random.random() < p_autosomal:
            female_T_alleles += 1

    male_freq_autosomal = male_T_alleles / (2 * n_males)
    female_freq_autosomal = female_T_alleles / (2 * n_females)

    print("--- Case 1: Autosomal Marker ---")
    print("This marker is present on a non-sex chromosome.")
    print(f"Frequency of allele 'T' in males: {male_freq_autosomal:.4f}")
    print(f"Frequency of allele 'T' in females: {female_freq_autosomal:.4f}")
    print("Result: Allele frequencies are very similar, leading to low Fst (close to 0).\n")


    # --- Case 2: Y-linked Marker ---
    # This marker is on the Y chromosome and thus exists ONLY in males.
    # Females do not have a Y chromosome.
    # This leads to maximum possible differentiation.
    p_y_linked = 0.6  # Frequency of allele 'M' on the Y chromosome
    
    male_M_alleles = 0
    
    # Simulate genotypes for males (1 allele each on the Y chromosome)
    for _ in range(n_males):
        if random.random() < p_y_linked:
            male_M_alleles += 1
    
    male_freq_y_linked = male_M_alleles / n_males
    # Females have no Y chromosome, so the frequency of any Y-linked allele is 0.
    female_freq_y_linked = 0.0

    print("--- Case 2: Y-linked Marker ---")
    print("This marker is present only on the Y chromosome.")
    print(f"Frequency of allele 'M' in males: {male_freq_y_linked:.4f}")
    print(f"Frequency of allele 'M' in females: {female_freq_y_linked:.4f}")
    print("Result: Allele frequencies are completely different, leading to high Fst (equal to 1).")
    print("\nThis simulation demonstrates how sex chromosome systems directly cause pronounced genetic differentiation for specific markers between males and females.")

# Run the simulation
simulate_population_genetics()