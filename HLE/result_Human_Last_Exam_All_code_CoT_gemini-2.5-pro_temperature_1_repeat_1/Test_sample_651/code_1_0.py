import sympy

def solve_billiards_limit():
    """
    Solves the billiards trajectory problem by finding the limit of M(theta) as theta -> 0.
    """
    # Step 1: Define theta as a symbol
    theta = sympy.Symbol('theta', real=True, positive=True)

    # Step 2: Formulate the maximization problem
    # M(theta) = sup(alpha). alpha is the angle of incidence.
    # alpha = pi/2 - beta, where beta is the angle between the trajectory and the side A.
    # Maximizing alpha is equivalent to minimizing beta.
    # M(theta) = pi/2 - beta_min
    # beta_min = arccos(sup(|cos(beta)|))
    # Therefore, M(theta) = arcsin(sup(|cos(beta)|))

    # Step 3: Analyze the boundary cases to find sup(|cos(beta)|)
    # The vector along side A is d_A = 5*exp(I*theta) - 5
    # The four extremal trajectories connect the endpoints of the unit arc {1, exp(I*theta)}
    # and the side A {5, 5*exp(I*theta)}.
    #
    # Through calculation, we find that the supremum of |cos(beta)| occurs for the
    # trajectories connecting (x=1, y=5e^(i*theta)) or (x=e^(i*theta), y=5).
    #
    # For these cases, |cos(beta)| is given by the expression:
    # sqrt(18*(1-cos(theta)) / (26-10*cos(theta)))

    sup_cos_beta = sympy.sqrt(18 * (1 - sympy.cos(theta)) / (26 - 10 * sympy.cos(theta)))

    # Step 4: Define M(theta)
    M_theta = sympy.asin(sup_cos_beta)

    print("The expression for the supremum angle of incidence M(theta) is:")
    print("M(theta) = arcsin(sqrt(18*(1 - cos(theta)) / (26 - 10*cos(theta))))")
    
    # As a check, let's find the behavior of M(theta) for small theta.
    # For small theta, sup(|cos(beta)|) is approximately 3*theta/4.
    # So, M(theta) is approximately arcsin(3*theta/4) ~= 3*theta/4.
    M_theta_div_theta = M_theta / theta
    limit_ratio = sympy.limit(M_theta_div_theta, theta, 0)
    print("\nThe ratio M(theta)/theta approaches a constant as theta -> 0:")
    print(f"lim_{{theta->0}} M(theta)/theta = {limit_ratio}")
    # The final equation is M(theta) / theta = 3/4
    print("So, for small theta, M(theta) is approximately (3/4) * theta")
    print("3 / 4 = 0.75")


    # Step 5: Compute the final limit of M(theta) as theta -> 0
    final_limit = sympy.limit(M_theta, theta, 0)

    print("\nThe limit of M(theta) as theta approaches 0 is:")
    print(f"lim_{{theta->0}} M(theta) = {final_limit}")

solve_billiards_limit()
<<<0>>>