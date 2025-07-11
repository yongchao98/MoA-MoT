import sympy

def print_overlap_integral_2s_2s():
    """
    Calculates and prints the analytical expression for the overlap integral
    of two 2s orbitals in a diatomic molecule like H2+.
    """
    # Define the symbolic variables
    zeta = sympy.Symbol('ζ')  # Effective nuclear charge
    R = sympy.Symbol('R')      # Internuclear distance

    # Define the dimensionless variable w
    w = (zeta * R) / 2

    # The analytical formula for the 2s-2s overlap integral, S(w), is:
    # S(w) = exp(-w) * [1 + w + (2/3)w^2 + (1/6)w^3 - (1/12)w^4]
    # We will construct this expression term by term for clarity.

    term0 = 1
    term1 = w
    term2 = sympy.Rational(2, 3) * w**2
    term3 = sympy.Rational(1, 6) * w**3
    term4 = -sympy.Rational(1, 12) * w**4

    # The polynomial part of the expression
    polynomial_part = term0 + term1 + term2 + term3 + term4

    # The full overlap integral expression
    S = sympy.exp(-w) * polynomial_part

    # Print the final equation clearly
    print("The overlap integral S for two 2s orbitals is given by:")
    print("S(ζ, R) = exp(-w) * P(w)")
    print("where w = ζ*R/2 and P(w) is a polynomial.")
    print("\nThe full expression is:")
    # Use pretty print for a nicely formatted equation
    final_equation_str = f"S(ζ, R) = exp(-({w})) * ({polynomial_part})"
    
    # To satisfy the "output each number" constraint, we can print the components.
    print(f"S(ζ, R) = exp(-ζ*R/2) * (c0 + c1*w + c2*w^2 + c3*w^3 + c4*w^4)")
    print(f"c0 = {term0}")
    print(f"c1 = {1}")
    print(f"c2 = {sympy.Rational(2, 3)}")
    print(f"c3 = {sympy.Rational(1, 6)}")
    print(f"c4 = {sympy.Rational(-1, 12)}")
    print("\nSubstituting w = ζ*R/2, the equation is:")

    # Print the equation with w substituted
    final_expr = sympy.exp(-w) * (1 + w + sympy.Rational(2, 3)*w**2 + sympy.Rational(1, 6)*w**3 - sympy.Rational(1, 12)*w**4)
    print(sympy.pretty(final_expr, use_unicode=True))


if __name__ == '__main__':
    print_overlap_integral_2s_2s()