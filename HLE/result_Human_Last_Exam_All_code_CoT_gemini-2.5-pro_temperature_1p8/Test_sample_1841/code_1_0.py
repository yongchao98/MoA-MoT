def calculate_and_explain_simpsons_index(sample_name, species_counts):
    """
    Calculates Simpson's Diversity Index and explains the steps.

    Args:
        sample_name (str): The name of the sample for printing.
        species_counts (list): A list where each number is the count of individuals
                               for a different species.
    """
    print(f"--- Analysis for {sample_name} ---")

    if not species_counts or sum(species_counts) < 2:
        print("Sample must contain at least two individuals to calculate the index.")
        print("-" * 35 + "\n")
        return

    # N = total number of organisms
    N = sum(species_counts)

    # sum_n_minus_1 = sum of n(n-1) for each species
    sum_n_n_minus_1_calc = [f"{n}*({n}-1)" for n in species_counts]
    sum_n_n_minus_1 = sum([n * (n - 1) for n in species_counts])

    # N(N-1)
    N_N_minus_1 = N * (N - 1)

    print(f"Species Counts (n): {species_counts}")
    print(f"Total Individuals (N): {N}")
    print(f"Sum of n(n-1) = {' + '.join(sum_n_n_minus_1_calc)} = {sum_n_n_minus_1}")
    print(f"N(N-1) = {N}*({N}-1) = {N_N_minus_1}")
    print("\nFinal Equation:")
    print(f"D = 1 - (Σn(n-1) / N(N-1))")
    print(f"D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")

    if N_N_minus_1 == 0:
        # This case is ill-defined but ecologically represents zero diversity.
        D = 0.0
        print("D = 0.0 (by definition, as there is no diversity)")
    else:
        fraction = sum_n_n_minus_1 / N_N_minus_1
        D = 1 - fraction
        print(f"D = 1 - {fraction:.4f}")
        print(f"Simpson's Diversity Index (D): {D:.4f}")
    
    print("-" * 35 + "\n")


# Case 1: A sample with multiple bat species (ecologically expected)
diverse_sample = [40, 25, 15, 10]
calculate_and_explain_simpsons_index("Diverse Sample (Multiple Species)", diverse_sample)

# Case 2: A sample where only one species was found (matches the student's result of D=0)
single_species_sample = [90]
calculate_and_explain_simpsons_index("Student's Sample (One Species)", single_species_sample)