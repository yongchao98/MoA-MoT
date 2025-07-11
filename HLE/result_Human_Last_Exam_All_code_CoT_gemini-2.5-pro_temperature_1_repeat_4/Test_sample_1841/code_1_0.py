def calculate_simpsons_diversity(species_counts):
    """
    Calculates Simpson's Diversity Index (D = 1 - sum[n(n-1)]/[N(N-1)])
    and explains the mathematical and ecological validity of a D=0 result.

    Args:
        species_counts (dict): A dictionary where keys are species names
                               and values are the number of individuals (n).
    """
    # Check if there are any individuals in the sample
    if not species_counts:
        print("The sample is empty. Cannot calculate diversity.")
        return

    # Get the list of individual counts (n for each species)
    n_values = list(species_counts.values())
    
    # Calculate N, the total number of individuals
    N = sum(n_values)

    print(f"Scenario: The student's sample contains only one species.")
    print(f"Sample data: {species_counts}")
    print("-" * 30)

    # Simpson's index is undefined for N < 2
    if N < 2:
        print(f"Total individuals (N) is {N}. Index is undefined for N < 2.")
        return

    # Calculate the numerator: sum of n(n-1) for all species
    numerator_sum = sum([n * (n - 1) for n in n_values])
    
    # Calculate the denominator: N(N-1)
    denominator = N * (N - 1)

    print("Calculating the Simpson's Diversity Index: D = 1 - (Σn(n-1) / N(N-1))")
    print("\n--- Equation components ---")
    
    # Print the breakdown of the calculation
    print(f"Total individuals (N) = {N}")
    for species, n in species_counts.items():
        print(f"Number of individuals for '{species}' (n) = {n}")
        print(f"Term n(n-1) for this species = {n} * ({n}-1) = {n * (n-1)}")
        
    print(f"\nSum of all n(n-1) terms (Σn(n-1)) = {numerator_sum}")
    print(f"Total term N(N-1) = {N} * ({N}-1) = {denominator}")
    
    # Avoid division by zero, although our N<2 check handles this
    if denominator == 0:
        diversity_index = 0 # Or undefined, but for a single species with N>=2, this won't be 0
    else:
        diversity_index = 1 - (numerator_sum / denominator)

    print("\n--- Final Calculation ---")
    print(f"D = 1 - ({numerator_sum} / {denominator})")
    print(f"D = 1 - {numerator_sum / denominator}")
    print(f"D = {diversity_index}")

    print("\n--- Conclusion ---")
    print("The calculation shows that a diversity index of 0 is mathematically possible. It occurs when a sample contains only one species.")
    print("However, since it is known that several bat species exist on the island, a result of D=0 indicates the sample was not representative of the true ecosystem. Therefore, the result is mathematically valid, but not ecologically valid.")


# --- Main execution ---
# Let's assume the student sampled for a month and found 50 bats, but they were all of the same species.
one_species_sample = {
    "Common Pipistrelle": 50
}

calculate_simpsons_diversity(one_species_sample)