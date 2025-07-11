def solve_phi4_fixed_point():
    """
    This script derives and explains the leading order expression for the 
    Wilson-Fisher fixed point coupling u* in phi^4 theory near d=4 dimensions.
    """
    print("This program calculates the leading order expression for the fixed point coupling u*")
    print("in the phi^4 theory in d = 4 - ε dimensions.")
    print("-" * 70)

    # Step 1: State the beta function for the coupling u.
    print("Step 1: The one-loop beta function, β(u), describes the change of the")
    print("coupling 'u' with the energy scale. It is given by:")
    print("\n  β(u) = -ε * u + C * u^2\n")
    print("where ε is a small parameter representing the deviation from 4 dimensions,")
    print("and C is a constant derived from loop calculations in the theory.")
    
    constant_numerator = 3.0
    constant_denominator = 16.0
    
    print(f"The value of C is ({int(constant_numerator)} / ({int(constant_denominator)} * π^2)).")
    print("-" * 70)

    # Step 2: Set the beta function to zero to find the fixed point.
    print("Step 2: A fixed point, u*, is a value of the coupling where the system")
    print("becomes scale-invariant. This occurs when the beta function is zero:")
    print("\n  β(u*) = 0\n")
    print("This gives the equation:")
    print(f"  -ε * u* + ({int(constant_numerator)} / ({int(constant_denominator)} * π^2)) * (u*)^2 = 0")
    print("-" * 70)

    # Step 3: Solve for the non-trivial (Wilson-Fisher) fixed point.
    print("Step 3: We are interested in the non-trivial fixed point where u* is not zero.")
    print("We can divide the equation by u* to get:")
    print(f"\n  -ε + ({int(constant_numerator)} / ({int(constant_denominator)} * π^2)) * u* = 0\n")
    print("Solving for u*, we find:")
    print("  u* = ε / C")
    print(f"  u* = ε / ({int(constant_numerator)} / ({int(constant_denominator)} * π^2))")
    print("-" * 70)

    # Step 4: Present the final expression.
    print("Step 4: The final leading order expression for the fixed point coupling u* is:")
    print("\n          (16 * π^2)")
    print(f"  u* = ————————— * ε\n")
    print("            3\n")

    print("To explicitly show each number in the final equation:")
    print(f"u* = ({int(constant_denominator)} * π^2 / {int(constant_numerator)}) * ε")


if __name__ == "__main__":
    solve_phi4_fixed_point()
