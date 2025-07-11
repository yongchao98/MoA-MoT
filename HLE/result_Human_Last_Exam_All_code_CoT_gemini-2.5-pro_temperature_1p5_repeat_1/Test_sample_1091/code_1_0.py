import sympy

def solve_limit_problem():
    """
    This script calculates the limit of n*P(n) as n -> infinity.
    It follows the plan outlined above, using symbolic mathematics to ensure precision.
    """
    # Step 1: Define problem parameters symbolically
    k = sympy.Symbol('k', positive=True)

    # There are 2k vectors of each type
    num_vectors_per_type = 2 * k

    # The three vector types
    vA = (sympy.S(1), sympy.S(0))
    vB = (sympy.S(1)/2, sympy.sqrt(3)/2)
    vC = (-sympy.S(1)/2, sympy.sqrt(3)/2)

    # Condition on the sum vector S: ||S||^2 <= R_squared
    R_squared = 2

    # Step 2: Calculate the components of the covariance matrix Sigma
    # Var(Sx) = sum of v_x^2 over all n vectors (since Var(epsilon_i)=1)
    var_sx = num_vectors_per_type * (vA[0]**2 + vB[0]**2 + vC[0]**2)
    var_sy = num_vectors_per_type * (vA[1]**2 + vB[1]**2 + vC[1]**2)
    cov_sxy = num_vectors_per_type * (vA[0]*vA[1] + vB[0]*vB[1] + vC[0]*vC[1])

    # Simplify the expressions
    var_sx_simp = sympy.simplify(var_sx)
    var_sy_simp = sympy.simplify(var_sy)
    cov_sxy_simp = sympy.simplify(cov_sxy)
    
    print("--- Step-by-Step Calculation ---")
    print("\n1. Calculate the Covariance Matrix of S=(Sx, Sy):")
    print(f"   Var(Sx) = 2k * (({vA[0]})**2 + ({vB[0]})**2 + ({vC[0]})**2) = {var_sx_simp}")
    print(f"   Var(Sy) = 2k * (({vA[1]})**2 + ({vB[1]})**2 + ({vC[1]})**2) = {var_sy_simp}")
    print(f"   Cov(Sx,Sy) = 2k * (({vA[0]}*{vA[1]}) + ({vB[0]}*{vB[1]}) + ({vC[0]}*{vC[1]})) = {cov_sxy_simp}")
    print(f"   So, the covariance matrix is Sigma = [[{var_sx_simp}, {cov_sxy_simp}], [{cov_sxy_simp}, {var_sy_simp}]]")

    # The resulting distribution is N(0, sigma^2 * I) where sigma^2 = 3k
    sigma_squared = var_sx_simp 

    # Step 3: Calculate the probability P(n) for large n.
    # For a bivariate normal N(0, sigma^2*I), P(||S||^2 <= R^2) = 1 - exp(-R^2 / (2*sigma^2))
    prob_in_k = 1 - sympy.exp(-R_squared / (2 * sigma_squared))
    
    print("\n2. Calculate the probability P(n):")
    print(f"   For large n, S follows a bivariate normal distribution. The probability P(||S||^2 <= {R_squared}) is given by:")
    print(f"   P(n) = 1 - exp(-R_squared / (2 * sigma^2)), where R_squared = {R_squared} and sigma^2 = {sigma_squared}")
    print(f"   P(n) = {prob_in_k}")

    # Step 4: Compute the limit of n*P(n)
    n_sym = sympy.Symbol('n', positive=True)
    # Substitute k = n/6 into the probability expression
    prob_in_n = prob_in_k.subs(k, n_sym/6)
    expression = n_sym * prob_in_n
    
    print("\n3. Calculate the limit of n*P(n):")
    print(f"   With n = 6k, we substitute k = n/6 into P(n). The expression to evaluate is:")
    print(f"   lim_{{n->oo}} ( n * P(n) ) = lim_{{n->oo}} ( {expression} )")
    
    # Calculate the limit as n -> infinity
    limit_val = sympy.limit(expression, n_sym, sympy.oo)
    
    print("\n--- Final Answer ---")
    print("The value of the limit is:")
    print(int(limit_val))

solve_limit_problem()
<<<2>>>