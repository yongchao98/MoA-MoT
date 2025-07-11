import sympy

def solve_billiards_limit():
    """
    Solves the billiards trajectory problem by finding the limit of M(theta) as theta -> 0.
    
    The solution involves analyzing the extremal trajectory that maximizes the angle alpha.
    This trajectory is found to be from x=1 to y=5*exp(i*theta).
    
    1. Define theta as a symbol.
    2. Define the expression for cos(alpha) for the extremal trajectory.
    3. Compute the Taylor series for cos(M(theta)) around theta=0.
    4. From the series, determine the behavior of M(theta) for small theta.
    5. Compute the limit of M(theta) as theta -> 0.
    """
    
    # Define theta as a symbolic variable
    theta = sympy.Symbol('theta', real=True, positive=True)

    # For the extremal trajectory from x=1 to y=5*exp(i*theta),
    # the cosine of the angle alpha=M(theta) can be expressed as:
    # cos(alpha) = (5*cos(3*theta/2) - cos(theta/2)) / sqrt(26 - 10*cos(theta))
    
    cos_M = (5 * sympy.cos(3 * theta / 2) - sympy.cos(theta / 2)) / sympy.sqrt(26 - 10 * sympy.cos(theta))

    # For small M, cos(M) is approximately 1 - M**2 / 2.
    # So, M**2 is approximately 2 * (1 - cos(M)).
    M_squared_approx = 2 * (1 - cos_M)

    # Compute the Taylor series expansion of M_squared_approx around theta = 0.
    # We are interested in the first non-zero term.
    series = sympy.series(M_squared_approx, theta, 0, 4)

    # The series is (49/16)*theta**2 + O(theta**4).
    # This means M(theta)^2 is approximately (49/16)*theta^2.
    # So, M(theta) is approximately sqrt(49/16)*theta = (7/4)*theta.
    
    m_theta_coeff_num = 7
    m_theta_coeff_den = 4
    
    # We found that M(theta) is proportional to theta for small theta.
    # M(theta) = (7/4) * theta
    
    # The final equation for M(theta) near theta=0 is:
    print(f"For small theta, M(theta) is approximately ({m_theta_coeff_num}/{m_theta_coeff_den}) * theta.")
    
    # Now we find the limit of M(theta) as theta goes to 0.
    limit_val = sympy.limit(m_theta_coeff_num * theta / m_theta_coeff_den, theta, 0)
    
    print(f"The limit of M(theta) as theta -> 0 is: {limit_val}")

solve_billiards_limit()