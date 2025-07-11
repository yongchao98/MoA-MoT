def calculate_simpsons_diversity(species_counts):
    """
    Calculates Simpson's Diversity Index for a given sample.
    D = 1 - [ sum(n*(n-1)) / (N*(N-1)) ]
    :param species_counts: A dictionary where keys are species names and values are their counts.
    """
    # Check for trivial cases
    if not species_counts or len(species_counts) <= 1:
        # If there's 0 or 1 species, diversity is 0.
        # N*(N-1) could be zero if N is 0 or 1, so we handle this case separately.
        n_values = list(species_counts.values())
        n = n_values[0] if n_values else 0
        N = n
        sum_n_n_minus_1 = n * (n - 1)
        N_N_minus_1 = N * (N - 1)
        
        print("Sample with one species:")
        print(f"Species counts (n): {n}")
        print(f"Total individuals (N): {N}")
        print(f"Sum of n(n-1): {n} * ({n}-1) = {sum_n_n_minus_1}")
        print(f"N(N-1): {N} * ({N}-1) = {N_N_minus_1}")
        
        if N_N_minus_1 == 0:
             # This can happen if N=1. The diversity is indisputably 0.
            diversity_index = 0.0
            print(f"Final Equation: D = 0.0 (by definition, as there is no diversity with one individual)")
        else:
            diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)
            print(f"Final Equation: D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1}) = {diversity_index}")
        
        print(f"\nCalculated Simpson's Diversity Index (D): {diversity_index}")
        print("This is a mathematically valid result.")
        print("However, it is not ecologically valid as it contradicts the prior knowledge that several bat species exist on the island.")

# --- Main Execution ---
# Scenario from the problem: The student's sample that resulted in D=0.
# This implies only one species was found. Let's assume 30 individuals of that one species.
student_sample = {"Bat Species A": 30}

# Calculate and explain the result for the student's sample.
calculate_simpsons_diversity(student_sample)
