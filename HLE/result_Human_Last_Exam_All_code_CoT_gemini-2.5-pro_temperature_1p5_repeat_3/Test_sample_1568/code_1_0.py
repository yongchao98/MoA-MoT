def solve_infinite_product():
    """
    This function formats and prints the closed-form expression for the given infinite product.
    The formula is derived using the Weierstrass product representation of the Gamma function.
    """
    # Define unicode characters for better mathematical notation
    PI = "\u03A0"
    GAMMA = "\u0393"
    OMEGA = "\u03C9"
    INFINITY = "\u221E"
    SUP_2 = "\u00B2"
    SUP_3 = "\u00B3"
    PI_SMALL = "\u03C0"
    SQRT = "\u221A"

    # Construct the left-hand side of the equation
    lhs = f"{PI}_{{n=3}}^{{{INFINITY}}} (1 - z{SUP_3}/n{SUP_3})"

    # Construct the right-hand side of the equation
    denominator = (f"(1 - z{SUP_3})(1 - z{SUP_3}/8)"
                   f"{GAMMA}(1 - z){GAMMA}(1 - {OMEGA}z){GAMMA}(1 - {OMEGA}{SUP_2}z)")
    rhs = f"1 / [{denominator}]"

    # Print the final equation
    print(f"{lhs} = {rhs}")
    print() # Add a newline for spacing

    # Print definitions of the special symbols used in the formula
    print("where:")
    print(f"  {GAMMA}(x) is the Gamma function.")
    print(f"  {OMEGA} = e^(2{PI_SMALL}i/3) = -1/2 + i*({SQRT}3)/2 is the principal cube root of unity.")

solve_infinite_product()