def calculate_simpsons_diversity(species_counts):
    """
    Calculates Simpson's Diversity Index and prints the steps.
    
    Args:
        species_counts (list): A list of integers where each integer is the number 
                               of individuals of a particular species.
    """
    N = sum(species_counts)
    
    if N < 2:
        print("Cannot calculate Simpson's Index for a total population (N) less than 2.")
        return

    sum_n_n_minus_1 = sum([n * (n - 1) for n in species_counts])
    N_N_minus_1 = N * (N - 1)
    
    if N_N_minus_1 == 0:
        # This case is unlikely if N>=2 but good practice to check
        simpson_index_D = 1.0 
    else:
        dominance_index = sum_n_n_minus_1 / N_N_minus_1
        simpson_index_D = 1 - dominance_index

    print(f"Sample counts (n for each species): {species_counts}")
    print(f"Total individuals (N): {N}")
    print(f"Sum of n(n-1) for each species (Σn(n-1)): {sum_n_n_minus_1}")
    print(f"Total individuals N(N-1): {N_N_minus_1}")
    print("\nCalculating the index: D = 1 - [Σn(n-1) / N(N-1)]")
    print(f"D = 1 - [{sum_n_n_minus_1} / {N_N_minus_1}]")
    print(f"D = 1 - {round(dominance_index, 4)}")
    print(f"Final Simpson's Diversity Index (D): {round(simpson_index_D, 4)}\n")


print("--- Case 1: The student's sample (hypothetical numbers) ---")
# This represents the student's sample, which must have had only one species to get D=0.
# Let's assume they found 30 bats, all of the same species.
student_sample = [30]
calculate_simpsons_diversity(student_sample)

print("--- Case 2: A more diverse sample (hypothetical numbers) ---")
# A hypothetical diverse sample with two species for comparison.
diverse_sample = [15, 15]
calculate_simpsons_diversity(diverse_sample)

print("--- Final Conclusion ---")
print("As demonstrated in Case 1, a Simpson's Diversity Index of 0 is a mathematically possible result.")
print("It occurs when the sample contains only one species (zero diversity).")
print("However, the problem states that multiple bat species are known to exist on the island.")
print("Therefore, the student's sample was not representative of the true ecological situation.")
print("The result is mathematically valid, but not ecologically valid.")
