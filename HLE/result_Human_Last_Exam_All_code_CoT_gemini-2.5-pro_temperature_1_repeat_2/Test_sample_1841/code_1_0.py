import math

def calculate_simpsons_index_from_counts(species_counts):
    """
    Calculates Simpson's Diversity Index from a dictionary of species counts.
    Prints the step-by-step calculation.
    """
    print("Analyzing the student's result using a hypothetical sample.")
    print("Let's assume the student sampled 42 bats, but they were all of a single species.\n")

    # n = number of individuals of a particular species
    # N = total number of individuals
    n_values = list(species_counts.values())
    N = sum(n_values)

    print(f"Species data from sample: {species_counts}")
    print(f"Total number of individuals (N) = {N}\n")

    # Calculate Σn(n-1)
    # This is the sum of n*(n-1) for each species
    sum_n_n_minus_1_terms = [n * (n - 1) for n in n_values]
    sum_n_n_minus_1 = sum(sum_n_n_minus_1_terms)

    print("Step 1: Calculate the numerator of the dominance fraction: sum[n(n-1)]")
    # Show the calculation for each term, even though there's only one
    for i, n in enumerate(n_values):
      print(f"  - Term for species {i+1}: {n} * ({n} - 1) = {sum_n_n_minus_1_terms[i]}")
    print(f"Numerator value = {sum_n_n_minus_1}\n")

    # Calculate N(N-1)
    if N < 2:
        print("Simpson's Index is not defined for N < 2.")
        return
        
    N_N_minus_1 = N * (N - 1)
    print("Step 2: Calculate the denominator of the dominance fraction: N(N-1)")
    print(f"  - Calculation: {N} * ({N} - 1) = {N_N_minus_1}\n")

    # Calculate dominance index
    dominance_index = sum_n_n_minus_1 / N_N_minus_1
    print("Step 3: Calculate the Simpson's Dominance Index (D_dominance)")
    print(f"  - D_dominance = Numerator / Denominator = {sum_n_n_minus_1} / {N_N_minus_1} = {dominance_index}\n")

    # Calculate diversity index (1 - D)
    diversity_index = 1 - dominance_index
    print("Step 4: Calculate the final Simpson's Diversity Index (1 - D_dominance)")
    print(f"  - Final Equation: 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")
    print(f"  - Calculation: 1 - {dominance_index} = {diversity_index}\n")

    print("Conclusion:")
    print("The calculation shows that a Simpson's Diversity Index of 0 is mathematically possible.")
    print("It occurs when the sample contains only one species.")
    print("However, this result contradicts the known ecological fact that the island has several bat species.")
    print("Therefore, the index value is mathematically valid for the sample collected, but it is not an ecologically valid representation of the entire island's diversity.")


# Hypothetical sample data where the student only found one species of bat.
# For example, 42 individuals of the species 'Myotis mystacinus'.
student_sample = {
    'Myotis mystacinus': 42
}

calculate_simpsons_index_from_counts(student_sample)
