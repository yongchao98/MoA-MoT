def calculate_simpsons_diversity(species_counts, scenario_name):
    """
    Calculates and explains Simpson's Diversity Index for a given sample.

    Args:
        species_counts: A list of counts for each species.
        scenario_name: A string describing the scenario.
    """
    print(f"--- {scenario_name} ---")
    
    if not species_counts or sum(species_counts) < 2:
        print("Sample size is too small to calculate diversity.\n")
        return

    # N = total number of all individuals
    N = sum(species_counts)
    
    # Σn(n-1) = sum of n*(n-1) for each species
    sum_n_n_minus_1 = sum([n * (n - 1) for n in species_counts])
    
    # N(N-1)
    N_N_minus_1 = N * (N - 1)

    # Simpson's Diversity Index D = 1 - [Σn(n-1) / N(N-1)]
    diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)

    print(f"Species counts: {species_counts}")
    print(f"Total individuals (N): {N}")
    print(f"Sum of n(n-1) for all species: {sum_n_n_minus_1}")
    print(f"Total individuals term N(N-1): {N_N_minus_1}")
    
    # Output the final equation with the numbers plugged in
    print("\nFinal Equation:")
    print(f"D = 1 - ( {sum_n_n_minus_1} / {N_N_minus_1} )")
    print(f"Calculated Simpson's Diversity Index: {diversity_index:.4f}\n")

# Scenario 1: Represents the student's sample data that would result in D=0
# This sample contains 60 individuals of a single bat species.
student_sample = [60]
calculate_simpsons_diversity(student_sample, "Student's Sample (One Species)")

# Scenario 2: A hypothetical sample reflecting the known diversity of the island
# This sample contains 30 individuals of two different bat species.
diverse_sample = [30, 30]
calculate_simpsons_diversity(diverse_sample, "Hypothetical Diverse Sample (Two Species)")