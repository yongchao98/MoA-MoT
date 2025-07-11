import math

def calculate_simpsons_diversity(species_counts, description):
    """
    Calculates the Simpson's Diversity Index and prints the step-by-step equation.
    
    Args:
    species_counts (list): A list of counts for each species found.
    description (str): A description of the scenario being calculated.
    """
    print(f"--- {description} ---")
    
    # Calculate N, the total number of individuals
    N = sum(species_counts)
    
    # The formula is undefined for N < 2.
    if N < 2:
        print("Simpson's Diversity Index is not applicable for a total sample size of less than 2.")
        return

    # Calculate the numerator: sum of n(n-1) for each species
    sum_n_n_minus_1 = sum([n * (n - 1) for n in species_counts])
    
    # Calculate the denominator: N(N-1)
    N_N_minus_1 = N * (N - 1)

    # Calculate Simpson's Diversity Index (D)
    # Avoid division by zero, although handled by the N < 2 check.
    if N_N_minus_1 == 0:
         diversity_index = 0 # This case implies N=0 or N=1
    else:
        diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)

    # Print the explanation and the final equation
    print("The Simpson's Diversity Index is calculated with the formula: D = 1 - [sum(n(n-1)) / N(N-1)]")
    print(f"Total individuals (N) = {N}")
    print(f"Numerator sum(n(n-1)) = {sum_n_n_minus_1}")
    print(f"Denominator N(N-1) = {N_N_minus_1}")
    print("\nFinal Equation:")
    print(f"D = 1 - ({sum_n_n_minus_1} / {N_N_minus_1})")
    print(f"D = 1 - {sum_n_n_minus_1 / N_N_minus_1 if N_N_minus_1 > 0 else 0}")
    print(f"D = {diversity_index:.4f}\n")

# Scenario from the problem: Student obtains a diversity index of 0.
# This happens when only one species is found. Let's assume the student found 45 bats of a single species.
student_sample = [45]
calculate_simpsons_diversity(student_sample, "Student's Sample (One Species)")

# For comparison, a hypothetical sample with multiple species.
# Let's say 20 bats of species A, 15 of species B, and 10 of species C were found.
diverse_sample = [20, 15, 10]
calculate_simpsons_diversity(diverse_sample, "Hypothetical Diverse Sample")