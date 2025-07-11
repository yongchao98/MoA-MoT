def calculate_simpsons_diversity():
    """
    This function demonstrates the calculation of Simpson's Diversity Index (D)
    for a sample that contains only one species, leading to a D value of 0.
    """
    # Imagine a sample where 45 bats were caught, but they all belong to a single species.
    # The list represents the number of individuals for each species found.
    # In this case, there's only one number because only one species was found.
    species_counts = [45]
    
    print("Scenario: A sample is collected with individuals from only one species.")
    print(f"Number of individuals for the single species (n): {species_counts[0]}\n")
    
    # N is the total number of all individuals
    N = sum(species_counts)
    
    # Σn(n-1) is the sum of n*(n-1) for each species
    # Since there's only one species, the sum is just n*(n-1) for that species.
    sum_n_n_minus_1 = 0
    for n in species_counts:
        # This part of the formula shows the dominance. If a community is dominated
        # by one species, this value will be high.
        sum_n_n_minus_1 += n * (n - 1)

    # N(N-1) is the denominator
    # It will only be undefined if N is 0 or 1. Let's assume more than 1 animal was sampled.
    if N <= 1:
        print("Simpson's Index is not defined for a total population (N) of 1 or less.")
        return

    N_N_minus_1 = N * (N - 1)
    
    # Simpson's Dominance Index (Lambda)
    lambda_index = sum_n_n_minus_1 / N_N_minus_1

    # Simpson's Diversity Index (D)
    D = 1 - lambda_index

    # Print the equation steps
    print("--- Simpson's Diversity Index Calculation ---")
    print("Formula: D = 1 - (Σn(n-1) / N(N-1))")
    print(f"Total individuals (N) = {N}")
    print(f"Sum of n(n-1) for all species [Σn(n-1)] = {sum_n_n_minus_1}")
    print(f"N(N-1) = {N}({N}-1) = {N_N_minus_1}")
    
    # Here we show the full equation with the numbers plugged in
    print("\nFinal Equation:")
    print(f"D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")
    print(f"D = 1 - {lambda_index}")
    print(f"D = {D}")
    
    print("\n--- Analysis ---")
    print("Mathematically, a Simpson's Diversity Index of 0 is valid. It occurs when the sample contains only one species.")
    print("Ecologically, however, this result contradicts the provided information that the island has historically been home to several bat species. A finding of D=0 represents the diversity *of the sample*, but it is not a valid representation of the *island's ecosystem diversity* given the prior knowledge. Therefore, the index value is mathematically valid but not ecologically valid in this context.")

calculate_simpsons_diversity()