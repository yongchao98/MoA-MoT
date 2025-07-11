import sys

def solve_critical_exponent():
    """
    Determines and explains the value of the critical exponent nu (ν) for a G4-theory.

    The G4-theoretical framework describes the Ising universality class. The value of the
    critical exponent ν, which characterizes the scaling of the correlation length,
    is dependent on the spatial dimension 'd'.
    """

    print("Analyzing the critical exponent ν for the G₄-theory (Ising universality class).")
    print("The value of ν depends on the spatial dimension 'd':")
    print("-" * 50)
    print("1. For d >= 4 (the upper critical dimension):")
    print("   Mean-field theory is exact, giving ν = 1/2.")
    print("\n2. For d = 3 (a physically relevant dimension):")
    print("   The value is known from high-precision numerical results, ν ≈ 0.63012.")
    print("\n3. For d = 2:")
    print("   The model is exactly solvable (Onsager's solution), yielding a precise, analytical value.")
    print("-" * 50)

    print("\nGiven the prompt's request for 'the precise value', the exact analytical result")
    print("from the 2-dimensional case is the most appropriate answer.")

    # Define the components for the final "equation" as requested
    dimension = 2
    exponent_symbol = "ν"
    equals_sign = "="
    exponent_value = 1

    # Print the final result, showing each component
    print(f"\nIn d={dimension} dimensions, the precise value of the critical exponent is:")
    sys.stdout.write(exponent_symbol + " ")
    sys.stdout.write(equals_sign + " ")
    sys.stdout.write(str(exponent_value) + "\n")


if __name__ == "__main__":
    solve_critical_exponent()
