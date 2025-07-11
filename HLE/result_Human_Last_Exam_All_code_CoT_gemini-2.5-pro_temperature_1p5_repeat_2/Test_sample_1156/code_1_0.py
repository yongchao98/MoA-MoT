import sympy as sp

def solve():
    """
    This function provides the derived analytical approximation for the
    normalised density of the invariant measure and prints its components.
    """
    # Define the symbol for the variable x
    x = sp.Symbol('x')

    # Define the constants and exponents in the equation
    numerator_const = 3
    denominator_const = 8
    exponent1_num = -1
    exponent1_den = 2
    exponent2_num = 1
    exponent2_den = 2

    # Construct the equation using the derived numbers
    rho = (sp.Rational(numerator_const, denominator_const) *
           (x**sp.Rational(exponent1_num, exponent1_den) +
            x**sp.Rational(exponent2_num, exponent2_den)))

    # Print the full equation
    print("The approximated normalised density rho(x) is:")
    sp.pprint(rho, use_unicode=False)
    print("\n")

    # Print each number in the final equation as requested
    print("The numbers in the final equation are:")
    print(f"Constant Numerator: {numerator_const}")
    print(f"Constant Denominator: {denominator_const}")
    print(f"Exponent 1 Numerator: {exponent1_num}")
    print(f"Exponent 1 Denominator: {exponent1_den}")
    print(f"Exponent 2 Numerator: {exponent2_num}")
    print(f"Exponent 2 Denominator: {exponent2_den}")

solve()