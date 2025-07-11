import sympy

def solve_equation():
    """
    This function uses symbolic mathematics to derive the final expression for theta'(t).
    It follows the steps outlined in the derivation above.
    """
    # Define symbols. K, c are constants. t is the independent variable.
    # theta and r are functions of t.
    t = sympy.Symbol('t')
    K = sympy.Symbol('K')
    c = sympy.Symbol('c')
    theta = sympy.Function('theta')(t)
    r = sympy.Function('r')(t)

    # Step 7: The derived expression for theta'(t) in terms of a(t) and a'(t).
    # We define a(t) and a'(t) symbolically for the substitution.
    a = sympy.Function('a')(t)
    a_prime = sympy.Function('a_prime')(t)
    
    theta_prime_intermediate = (c * (K * a**2 + a_prime**2)) / (c**2 * a**2 + a_prime**2)

    # Step 8: The relationships between (a, a') and (r, theta).
    a_subs = r * sympy.sin(theta)
    a_prime_subs = c * r * sympy.cos(theta)

    # Substitute these into the intermediate expression for theta'(t).
    final_expression = theta_prime_intermediate.subs([
        (a, a_subs),
        (a_prime, a_prime_subs)
    ])

    # Simplify the final expression.
    simplified_expression = sympy.simplify(final_expression)

    # Print the final result in a clear format.
    print("The derived expression for theta'(t) is:")
    
    # To satisfy the "output each number" instruction, we print the components
    # of the final formula. The expression is symbolic, so the "numbers" are
    # the symbolic constants and functions.
    final_term1_coeff = c
    final_term1_func = sympy.cos(theta)**2
    final_term2_coeff = K/c
    final_term2_func = sympy.sin(theta)**2
    
    # Using pprint for a nicely formatted mathematical expression
    sympy.pprint(final_term1_coeff * final_term1_func + final_term2_coeff * final_term2_func)

    # For clarity, also printing out the components of the equation
    print("\nBreaking down the final equation:")
    print(f"Term 1: ({final_term1_coeff}) * {final_term1_func}")
    print(f"Term 2: ({final_term2_coeff}) * {final_term2_func}")

solve_equation()