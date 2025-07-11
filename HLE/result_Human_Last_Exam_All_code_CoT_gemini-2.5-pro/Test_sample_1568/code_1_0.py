def solve_infinite_product():
    """
    This function prints the closed-form expression for the infinite product:
    P = product_{n=3 to infinity} (1 - z^3 / n^3)
    
    The final equation is derived using the Weierstrass product representation
    of the Gamma function. The resulting formula involves Gamma functions and
    the complex cube roots of unity.
    """

    # The equation for the product from n=3 to infinity.
    # The left-hand side (LHS) of the equation.
    lhs = "\u220F_{n=3 to \u221E} (1 - z^3/n^3)"

    # The right-hand side (RHS) of the equation.
    # The numerator is 1.
    numerator = "1"
    
    # The denominator consists of the terms for n=1, n=2, and the product of Gamma functions.
    # We will first show the unsimplified terms with 1^3 and 2^3.
    denominator_unsimplified = "(1 - z^3/1^3) * (1 - z^3/2^3) * \u0393(1-z) * \u0393(1-z*e^(2\u03C0i/3)) * \u0393(1-z*e^(-2\u03C0i/3))"
    
    # And the simplified version where 1^3=1 and 2^3=8.
    denominator_simplified = "  (1 - z^3) * (1 - z^3/8) * \u0393(1-z) * \u0393(1-z*e^(2\u03C0i/3)) * \u0393(1-z*e^(-2\u03C0i/3))"

    # Print the equation in a structured way.
    print("The closed-form expression for the infinite product is:")
    print()
    print(f"{lhs} =")
    print()
    print(f" {numerator.center(80)}")
    print(" " + "-" * 80)
    print(f" {denominator_unsimplified}")
    print()
    print("After simplifying the terms for n=1 and n=2 (1^3=1, 2^3=8), the final expression is:")
    print()
    print(f"{lhs} =")
    print()
    print(f" {numerator.center(80)}")
    print(" " + "-" * 80)
    print(f"{denominator_simplified}")
    print()
    print("where \u0393 is the Gamma function, and i is the imaginary unit.")

solve_infinite_product()