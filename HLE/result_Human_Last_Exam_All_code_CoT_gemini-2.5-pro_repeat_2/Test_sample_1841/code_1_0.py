import math

def calculate_simpsons_diversity(species_counts):
    """
    Calculates Simpson's Diversity Index and prints the step-by-step calculation.

    Args:
        species_counts (dict): A dictionary where keys are species names and
                               values are the number of individuals of that species.
    """
    # The student's sample resulted in an index of 0.
    # This is only mathematically possible if only one species was found.
    # Let's create a hypothetical sample that reflects this.
    print("Hypothetical Sample Data based on student's result:")
    for species, n in species_counts.items():
        print(f"- {species}: {n} individuals (n)")

    # N = total number of organisms of all species
    N = sum(species_counts.values())

    # sum_n_minus_1 = sum of n(n-1) for all species
    sum_n_minus_1 = sum(n * (n - 1) for n in species_counts.values())

    # The denominator in the formula is N(N-1)
    # Handle the case where N < 2, which would make the index undefined or 0.
    if N < 2:
        dominance_index = 1
        diversity_index = 0
        denominator = "N/A (N < 2)"
    else:
        denominator = N * (N - 1)
        dominance_index = sum_n_minus_1 / denominator
        diversity_index = 1 - dominance_index

    print("\n--- Calculation Steps ---")
    print(f"The formula for Simpson's Diversity Index is: D = 1 - [ Σn(n-1) / N(N-1) ]")
    print(f"1. Total number of individuals (N) = {N}")

    # Show calculation for each n(n-1)
    n_calcs = [f"{n}({n}-1)" for n in species_counts.values()]
    n_results = [n * (n - 1) for n in species_counts.values()]
    print(f"2. Sum of n(n-1) for each species [Σn(n-1)] = {' + '.join(n_calcs)} = {sum(n_results)}")

    if N >= 2:
      print(f"3. Denominator N(N-1) = {N}({N}-1) = {denominator}")
    else:
      print(f"3. Denominator N(N-1) = {denominator}")


    print("\n--- Final Equation ---")
    if N >= 2:
      final_equation = f"D = 1 - ( {sum_n_minus_1} / {denominator} )"
      print(final_equation)
      print(f"D = 1 - {dominance_index}")
    else:
      print("D = 0 (as diversity is not possible with less than 2 individuals)")


    print("\n--- Result ---")
    print(f"Calculated Simpson's Diversity Index (D): {diversity_index:.4f}")

    print("\n--- Analysis ---")
    print("The result D=0 is mathematically valid, as shown by the calculation for a sample containing only one species.")
    print("However, it is not ecologically valid because we know from the problem description that multiple bat species exist on the island.")
    print("This means the student's sample was not representative of the true diversity of the entire island.")
    print("Therefore, the correct choice is 'Mathematically valid, but not ecologically valid.'")

# This is a hypothetical sample that would produce a diversity index of 0.
# It assumes the student found 30 bats, all of the same species.
student_sample = {"Common Noctule Bat": 30}

calculate_simpsons_diversity(student_sample)
