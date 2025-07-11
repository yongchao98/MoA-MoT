def simpsons_diversity_analysis():
    """
    Calculates and explains a scenario where Simpson's Diversity Index is 0.
    This demonstrates the mathematical validity but ecological invalidity
    described in the problem.
    """
    # Let's assume the student's sample consisted of 75 bats, all of the same species.
    # This represents a community with no diversity.
    # The list contains the count for each species found.
    species_counts = [75]

    # N is the total number of individuals from all species.
    N = sum(species_counts)

    # Σn(n-1) is the sum of n*(n-1) for each species.
    sum_n_n_minus_1 = 0
    for n in species_counts:
        sum_n_n_minus_1 += n * (n - 1)

    # N*(N-1) is the denominator part.
    N_N_minus_1 = N * (N - 1)

    print("--- Analysis of the Student's Result (D=0) ---")
    print(f"Scenario: A sample is collected with {N} individuals from only 1 species.")
    print(f"Species counts (n): {species_counts}")
    print(f"Total individuals (N): {N}\n")

    print("Formula: D = 1 - [Σn(n-1)] / [N(N-1)]\n")

    print("Step 1: Calculate the numerator Σn(n-1)")
    # We show the calculation for each species, even though there's only one.
    numerator_calc_str = " + ".join([f"{n}({n}-1)" for n in species_counts])
    print(f"Σn(n-1) = {numerator_calc_str} = {sum_n_n_minus_1}\n")

    print("Step 2: Calculate the denominator N(N-1)")
    print(f"N(N-1) = {N}({N}-1) = {N_N_minus_1}\n")

    print("Step 3: Calculate the final index D")
    if N_N_minus_1 > 0:
        diversity_index = 1 - (sum_n_n_minus_1 / N_N_minus_1)
        # Final equation with all numbers plugged in
        print(f"D = 1 - [{sum_n_n_minus_1}] / [{N_N_minus_1}]")
        print(f"D = 1 - {sum_n_n_minus_1 / N_N_minus_1}")
        print(f"D = {diversity_index}")
    else:
        # This case handles N=0 or N=1, where diversity is 0.
        diversity_index = 0
        print("D = 0 (Since total individuals N is less than 2, diversity is zero)")

    print("\n--- Conclusion ---")
    print("The calculation shows that D=0 is a mathematically possible result, occurring when the sample has zero diversity (only one species).")
    print("However, this contradicts the known fact that several bat species live on the island. Thus, the sample is not representative of the ecosystem, making the result ecologically invalid.")
    print("\nCorrect Answer: Mathematically valid, but not ecologically valid.")

if __name__ == "__main__":
    simpsons_diversity_analysis()