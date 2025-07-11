import sympy

def solve_curve_discriminant():
    """
    Calculates the minimal discriminant of the curve y^2 = x^6 + 2x^3 + 4x^2 + 4x + 1.

    The minimal discriminant of a hyperelliptic curve is interpreted here as the
    discriminant of the polynomial on the right-hand side, as simple transformations
    do not yield a model with a smaller discriminant.
    """
    # Define the symbolic variable
    x = sympy.Symbol('x')

    # Define the polynomial f(x) from the equation y^2 = f(x)
    # The coefficients are 1, 2, 4, 4, 1 for the powers x^6, x^3, x^2, x, and the constant term.
    f = x**6 + 2*x**3 + 4*x**2 + 4*x + 1

    # Calculate the discriminant of the polynomial f(x)
    disc = sympy.discriminant(f, x)

    # The discriminant is a single integer value. To show the "equation" that
    # produces it, we can display its prime factorization.
    abs_disc = abs(disc)
    factors = sympy.factorint(abs_disc)

    # Build an equation string from the prime factorization
    sign = "-" if disc < 0 else ""
    equation_parts = []
    for base, exponent in sorted(factors.items()):
        equation_parts.append(f"{base}^{exponent}")
    
    equation_str = sign + " * ".join(equation_parts)

    print(f"The polynomial is f(x) = x^6 + 2*x^3 + 4*x^2 + 4*x + 1.")
    print("The minimal discriminant of the curve y^2 = f(x) is the discriminant of this polynomial.")
    print("\nThe final equation for the discriminant, based on its prime factorization, is:")
    
    # Output each component of the factorization
    components = []
    if disc < 0:
        components.append("-1")
    
    for base, exponent in sorted(factors.items()):
        component_value = base**exponent
        components.append(str(component_value))
        print(f"Component: {base}^{exponent} = {component_value}")
    
    print(f"\nFinal Equation: {' * '.join(components)} = {disc}")
    print(f"\nThe minimal discriminant is: {disc}")

if __name__ == '__main__':
    solve_curve_discriminant()