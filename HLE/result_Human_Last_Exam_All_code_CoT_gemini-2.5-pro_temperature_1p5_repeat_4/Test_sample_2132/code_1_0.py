import math

def calculate_bottcher_complexity():
    """
    Calculates the Böttcher Molecular Complexity for the product of the
    Favorskii rearrangement of 2-chlorocyclohexanone.
    """
    # The product is cyclopentanecarboxylic acid (C6H10O2).
    # We first determine the parameters N, E, and P for the BMC formula.

    # N: Number of non-hydrogen atoms (6 Carbon + 2 Oxygen)
    N = 8

    # E: Number of bonds between non-hydrogen atoms.
    # (5 in the ring + 1 connecting the COOH group + 2 in the COOH group)
    E = 8

    # P: Number of partitions = C (components) + L (cycles)
    # C = 1 (the molecule is one connected component)
    # L = E - N + C = 8 - 8 + 1 = 1 (one ring)
    P = 1 + 1

    # --- Calculation ---
    # Calculate P factorial
    p_factorial = math.factorial(P)

    # Calculate the exponent term (N - E + P)
    exponent = N - E + P

    # Calculate the base term (E / N)
    base = E / N

    # Calculate the final BMC value
    bmc_value = N * (base ** exponent) / p_factorial

    # --- Output the results ---
    print("The product of the Favorskii rearrangement of 2-chlorocyclohexanone is cyclopentanecarboxylic acid.")
    print("The Böttcher Molecular Complexity (BMC) is calculated as: BMC = N * (E / N)^(N - E + P) / P!")
    print("\nFor cyclopentanecarboxylic acid:")
    print(f"N (non-hydrogen atoms) = {N}")
    print(f"E (bonds between non-hydrogen atoms) = {E}")
    print(f"P (partitions = components + cycles) = {P}")

    print("\nPlugging the values into the formula:")
    # Print the equation with all the numbers filled in
    print(f"BMC = {N} * ({E} / {N})^({N} - {E} + {P}) / {P}!")
    print(f"BMC = {N} * ({int(base)})^{exponent} / {p_factorial}")
    print(f"BMC = {N * (base ** exponent)} / {p_factorial}")
    print(f"\nFinal Böttcher Molecular Complexity = {bmc_value}")

# Execute the function
calculate_bottcher_complexity()