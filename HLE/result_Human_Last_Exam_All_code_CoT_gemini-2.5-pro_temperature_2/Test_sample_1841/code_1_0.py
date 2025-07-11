def calculate_simpsons_index(species_counts):
    """
    Calculates Simpson's Diversity Index (D).
    D = 1 - (sum(n(n-1)) / N(N-1))
    where:
    n = the number of individuals of each species
    N = the total number of individuals of all species
    """
    # Let's model the student's sample.
    # A result of D=0 implies only one species was found.
    # We'll create a hypothetical sample of 25 bats, all of a single species.
    # The dictionary represents {'species_name': count}
    if not species_counts:
        print("Sample is empty. Cannot calculate diversity.")
        return

    n_values = list(species_counts.values())
    N = sum(n_values)

    if N < 2:
        print(f"Total number of individuals (N={N}) is less than 2. Simpson's index is undefined.")
        # D would be 1 - (0/0) which is undefined. Ecologically, diversity requires at least 2 individuals.
        # But for D=0, N must be > 1.
        return

    # Calculate the sum of n(n-1) for each species
    sum_n_minus_1 = sum([n * (n - 1) for n in n_values])
    
    # Calculate N(N-1)
    N_N_minus_1 = N * (N - 1)

    # If the denominator is zero (i.e., N=0 or N=1), diversity is 0 or undefined.
    # Our previous check handles this, but it's good practice.
    if N_N_minus_1 == 0:
        dominance_index = 0
    else:
        dominance_index = sum_n_minus_1 / N_N_minus_1
        
    diversity_index = 1 - dominance_index

    # Print the explanation and the final equation
    print("Scenario: A sample contains individuals of only one species.")
    print(f"Sampled species counts: {species_counts}")
    print(f"The number of individuals of the single species (n) is {n_values[0]}.")
    print(f"The total number of individuals (N) is {N}.")
    print("\nCalculating the parts of the formula D = 1 - (sum(n(n-1)) / N(N-1)):")
    print(f"Numerator sum(n(n-1)): {n_values[0]} * ({n_values[0]} - 1) = {sum_n_minus_1}")
    print(f"Denominator N(N-1): {N} * ({N} - 1) = {N_N_minus_1}")
    
    print("\nFinal Equation:")
    print(f"D = 1 - ({sum_n_minus_1} / {N_N_minus_1})")
    print(f"D = 1 - {dominance_index}")
    print(f"D = {diversity_index}")
    
    print("\nConclusion: A Simpson's Diversity Index of 0 is mathematically valid and occurs when the sample contains only one species.")

# Run the calculation for a single-species sample.
# For example, 25 individuals of 'Common Pipistrelle'.
single_species_sample = {'Common Pipistrelle': 25}
calculate_simpsons_index(single_species_sample)
