def solve_inductance():
    """
    This script presents the derivation and final expression for the change
    in mutual inductance as described in the problem.
    """

    # The derivation outlined in the plan leads to the conclusion that the
    # change in mutual inductance ΔM is equal to the initial mutual inductance M1,
    # as the shells act as perfect shields making M2 = 0.

    # The calculation of M1 for two parallel wire pairs in the limit d >> h gives
    # the following result for the mutual inductance per unit length.

    print("The expression for the change in mutual inductance per unit length, ΔM, is:")
    final_expression = "(μ0 * h**2) / (2 * π * d**2)"
    print(final_expression)

    print("\nIn the final equation, the components are:")

    # To satisfy the prompt "output each number in the final equation",
    # we identify and print each symbolic and numeric component.
    
    term_1 = "μ0 (permeability of free space)"
    term_2 = "h (separation of wires in a circuit)"
    exponent_1 = "2"
    term_3 = "2 (a constant)"
    term_4 = "π (the constant pi)"
    term_5 = "d (separation between circuits)"
    exponent_2 = "2"

    print(f"Numerator part 1: {term_1}")
    print(f"Numerator part 2: {term_2}")
    print(f"Numerator part 2 is raised to the power of: {exponent_1}")
    
    print(f"Denominator part 1: {term_3}")
    print(f"Denominator part 2: {term_4}")
    print(f"Denominator part 3: {term_5}")
    print(f"Denominator part 3 is raised to the power of: {exponent_2}")

solve_inductance()