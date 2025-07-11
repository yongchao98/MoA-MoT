import sympy

def solve_sugra_parameters():
    """
    This function calculates the parameters alpha^2 and beta based on the
    provided supergravity problem.
    """

    # Part 1: Determine beta
    # From the cancellation of S-independent terms in the variation of L_cos,
    # we derived the equation: 1/2 - 2*beta = 0.
    beta = sympy.Symbol('beta')
    beta_eq = sympy.Eq(1/2 - 2 * beta, 0)
    beta_sol = sympy.solve(beta_eq, beta)[0]

    # Part 2: Determine alpha^2
    # From the bosonic potential and Einstein's equations for the AdS vacuum,
    # we derived the relation R = -3 * kappa^2 * alpha^2.
    R = sympy.Symbol('R')
    kappa = sympy.Symbol('kappa')
    alpha_sq = sympy.Symbol('alpha**2')
    
    # The equation is R = -3 * kappa^2 * alpha^2
    # We solve for alpha^2 in terms of R and kappa.
    alpha_sq_expr = -R / (3 * kappa**2)
    
    print("The value of beta is determined by requiring the S-independent part of the variation of L_cos to vanish.")
    print("The equation is: 1/2 - 2*beta = 0")
    print(f"Solving for beta, we find:")
    print(f"beta = {beta_sol}")
    
    print("\nThe value of alpha^2 is determined from the vacuum solution.")
    print("The scalar potential for S is minimized, leading to an AdS spacetime.")
    print("The relation from Einstein's equations is: R = -3 * kappa**2 * alpha**2")
    print(f"Solving for alpha**2, we get:")
    # We construct the equation string for alpha^2 manually for clear output
    equation_str = f"alpha**2 = {sympy.pretty(alpha_sq_expr, use_unicode=False)}"
    print(equation_str)
    
    # For the final answer format
    print("\n<<<")
    print(f"alpha**2 = {-R/(3*kappa**2)}")
    print(f"beta = {beta_sol}")
    print(">>>")


solve_sugra_parameters()