def calculate_and_explain_simpsons_index(sample_data, sample_name):
    """
    Calculates Simpson's Diversity Index and prints a step-by-step explanation.
    
    Args:
        sample_data (dict): A dictionary where keys are species names and values are their counts.
        sample_name (str): A name for the sample being analyzed.
    """
    print(f"--- Analyzing Sample: '{sample_name}' ---")
    
    counts = list(sample_data.values())
    
    # Calculate N, the total number of individuals
    N = sum(counts)
    if N <= 1:
        print("Simpson's Index cannot be calculated for a sample with 1 or 0 individuals.\n")
        return

    # Calculate Σn(n-1)
    sum_n_n_minus_1_terms = [n * (n - 1) for n in counts]
    sum_n_n_minus_1 = sum(sum_n_n_minus_1_terms)

    # Calculate N(N-1)
    N_N_minus_1 = N * (N - 1)
    
    # Calculate the Index D
    dominance_index = sum_n_n_minus_1 / N_N_minus_1
    diversity_index = 1 - dominance_index

    # Print the step-by-step calculation
    print("Formula: D = 1 - [Σn(n-1) / (N(N-1))]")
    print(f"Species counts (n): {counts}")
    print(f"Total individuals (N): {N}")
    print(f"Calculating the numerator Σn(n-1):")
    n_calculation_str = ' + '.join([f'{n}({n-1})' for n in counts])
    n_result_str = ' + '.join(map(str, sum_n_n_minus_1_terms))
    print(f"  Σn(n-1) = {n_calculation_str} = {n_result_str} = {sum_n_n_minus_1}")
    print(f"Calculating the denominator N(N-1):")
    print(f"  N(N-1) = {N}({N-1}) = {N_N_minus_1}")
    print(f"Final Equation:")
    print(f"  D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")
    print(f"  D = 1 - {dominance_index:.4f}")
    print(f"  D = {diversity_index:.4f}")
    print(f"Result for '{sample_name}': {diversity_index:.4f}\n")


# Scenario 1: A hypothetical sample representing known island diversity
multi_species_sample = {
    'Fruit Bat': 15,
    'Free-tailed Bat': 10,
    'Leaf-nosed Bat': 8
}
calculate_and_explain_simpsons_index(multi_species_sample, "Diverse Sample")

# Scenario 2: The student's sample that yielded an index of 0
# This requires only one species to be present.
student_sample = {
    'Fruit Bat': 33 # The total N from the other sample
}
calculate_and_explain_simpsons_index(student_sample, "Student's Actual Sample")

print("Conclusion:")
print("The code demonstrates that a Simpson's Diversity Index of 0 is a mathematically possible result, occurring when the sample contains only one species.")
print("However, this result is not ecologically valid for the island as a whole, because it is known that multiple bat species are present.")
print("This points to a discrepancy between the sample data and the true ecosystem diversity.")
