import sympy

def solve():
    """
    This function calculates the limit of M(theta) as theta approaches 0.
    """
    # Define the symbolic variable theta
    theta = sympy.symbols('theta')

    # As derived in the explanation, the cosine of the maximum angle M(theta) is given by:
    # cos(M(theta)) = 4*cos(theta/2) / sqrt(26 - 10*cos(theta))
    cos_M_theta = 4 * sympy.cos(theta/2) / sympy.sqrt(26 - 10 * sympy.cos(theta))

    # M(theta) is the arccos of this expression
    M_theta = sympy.acos(cos_M_theta)

    # We want to find the limit of M(theta) as theta approaches 0.
    # We can do this by first finding the limit of its argument.
    
    # Numerator limit
    num = 4 * sympy.cos(theta/2)
    limit_num = sympy.limit(num, theta, 0)
    
    # Denominator limit
    den = sympy.sqrt(26 - 10 * sympy.cos(theta))
    limit_den = sympy.limit(den, theta, 0)

    # Limit of the argument
    limit_arg = limit_num / limit_den

    # The final limit is the arccos of the argument's limit
    final_limit = sympy.acos(limit_arg)

    # Output the steps of the final calculation as requested
    print(f"To find the limit of M(theta), we first evaluate the limit of its cosine:")
    print(f"lim_{{theta->0}} cos(M(theta)) = (lim_{{theta->0}} 4*cos(theta/2)) / (lim_{{theta->0}} sqrt(26 - 10*cos(theta)))")
    print(f"The limit of the numerator is: 4 * cos(0) = {limit_num}")
    print(f"The limit of the denominator is: sqrt(26 - 10*cos(0)) = sqrt(26 - 10) = sqrt(16) = {limit_den}")
    print(f"So, the limit of cos(M(theta)) is: {limit_num} / {limit_den} = {limit_arg}")
    print(f"The limit of M(theta) is arccos of this value:")
    print(f"lim_{{theta->0}} M(theta) = arccos({limit_arg}) = {final_limit}")

solve()