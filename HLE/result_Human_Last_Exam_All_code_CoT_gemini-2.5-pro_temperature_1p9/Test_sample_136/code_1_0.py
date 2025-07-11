import sympy

def solve_probability_limit():
    """
    Calculates the limit of the conditional probability using a Bivariate Poisson model
    for the random walk local times, with parameters derived from Green's function asymptotics.
    """
    # Define symbols
    n = sympy.Symbol('n', real=True, positive=True)
    ln = sympy.log
    
    # Identify the vertex x0.
    # The origin is 0=(0,0). Its neighbours are (1,0), (-1,0), (0,1), (0,-1).
    # A vertex x0 with two common neighbours with 0 is, for instance, (1,1).
    # The common neighbours are (1,0) and (0,1).
    # The squared Euclidean distance ||x0||^2 is 1^2 + 1^2 = 2.
    norm_x0_sq = 2
    
    # Step 1: Define mean occupation time
    t_n = n**2 * ln(n)**2
    N = n**2
    mu_n = t_n / N
    print(f"Mean occupation time for any site, E[L(v)] = mu_n = {mu_n}")
    
    # Step 2: Asymptotic formulas for the Green's function G(x,y) on the torus.
    # The Green's function G(x,y) for a simple random walk on the torus Tn is asymptotically
    # G(x,y) ~ (1/pi) * (ln(n) - ln(||x-y||)).
    # We only need the ratio of Green's functions, so the (1/pi) constant cancels out.
    G_00_expr = ln(n) # Proportional to G(0,0)
    G_0x0_expr = ln(n) - ln(sympy.sqrt(norm_x0_sq)) # Proportional to G(0, x0)
    
    print("\nAsymptotic formulas for Green's function (up to a constant factor):")
    print(f"G(0,0) ~ {G_00_expr}")
    print(f"G(0,x0) ~ {G_0x0_expr}")

    # Step 3: Compute the correlation coefficient rho
    rho_n = G_0x0_expr / G_00_expr
    rho_n_simplified = sympy.simplify(rho_n)
    print(f"\nCorrelation between local times at 0 and x0, rho_n = {rho_n_simplified}")
    
    # Step 4: Model with Bivariate Poisson Distribution
    # L0 = Y0 + Y12, Lx0 = Yx0 + Y12
    # E[L0] = lambda_0 + lambda_12 = mu_n
    # Cov(L0, Lx0) = lambda_12
    # The model implies: lambda_12 = rho_n * mu_n
    # lambda_x0 = mu_n - lambda_12 = mu_n * (1 - rho_n)
    
    lambda_x0 = mu_n * (1 - rho_n)
    lambda_x0_simplified = sympy.simplify(lambda_x0)
    print(f"\nThe parameter lambda_x0 for the independent part of the Poisson model is:")
    print(f"lambda_x0 = mu_n * (1 - rho_n) = {lambda_x0_simplified}")
    
    # Step 5: Compute the conditional probability
    # P(Lx0=0 | L0=0) = P(Yx0+Y12=0 | Y0+Y12=0) = P(Yx0=0) = exp(-lambda_x0)
    prob_expr = sympy.exp(-lambda_x0)
    prob_expr_simplified = sympy.simplify(prob_expr)
    
    print(f"\nThe conditional probability P(x0 not visited | 0 not visited) is e^(-lambda_x0):")
    # To show each number, let's substitute the terms in the final expression
    base = n
    exponent = -ln(2)/2
    final_form = base**exponent
    print(f"P(Lx0=0|L0=0) = {prob_expr_simplified}")
    print(f"This expression can be written as: {base}^({exponent})")


    # Step 6: Compute the limit as n -> infinity
    final_limit = sympy.limit(prob_expr, n, sympy.oo)
    
    print("\nThe limit as n -> infinity is:")
    print(final_limit)

if __name__ == '__main__':
    solve_probability_limit()