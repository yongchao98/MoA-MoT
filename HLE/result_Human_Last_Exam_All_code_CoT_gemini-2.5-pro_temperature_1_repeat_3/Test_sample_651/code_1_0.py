import sympy

def solve_limit():
    """
    This function calculates the limit of M(theta) as theta -> 0 based on the
    derived expression M(theta) = 3*theta/4.
    """
    # Define theta as a symbol
    theta = sympy.Symbol('theta')

    # As derived in the thinking steps, for small positive theta,
    # the supremum of the angle alpha is given by the function M(theta).
    M_theta = 3 * theta / 4

    # We need to find the limit of M(theta) as theta approaches 0.
    limit_M = sympy.limit(M_theta, theta, 0)

    # The final equation is: lim_{theta->0} (3 * theta / 4) = 0
    # The problem asks to output each number in the final equation.
    numerator = 3
    denominator = 4
    
    print("The derived expression for M(theta) for small theta is:")
    # Using pprint for a clean mathematical display
    sympy.pprint(M_theta, use_unicode=True)
    
    print(f"\nThe limit is calculated as: lim_{{theta->0}} ({numerator}*theta / {denominator})")
    
    print(f"\nThe final result of the limit is: {limit_M}")

solve_limit()