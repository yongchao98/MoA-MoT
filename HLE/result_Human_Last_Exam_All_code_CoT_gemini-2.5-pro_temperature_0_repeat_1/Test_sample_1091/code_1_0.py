import sympy as sp

def solve_limit_problem():
    """
    This function solves the limit problem analytically using the Central Limit Theorem
    and symbolic computation with sympy.
    """
    # Step 1: Define symbolic variables for the problem.
    # n is the total number of vectors, and n = 6k.
    n = sp.Symbol('n', positive=True, real=True)
    k = sp.Symbol('k', positive=True, real=True)

    print("Step 1: Applying the Central Limit Theorem.")
    print("The sum S=(Sx, Sy) follows a bivariate normal distribution for large n.")
    print("The mean is (0, 0). We need to compute the covariance matrix.")

    # There are 2k vectors of each type: (1,0), (0.5, sqrt(3)/2), (-0.5, sqrt(3)/2).
    # Let X_A, X_B, X_C be the sums of Rademacher variables for each vector type.
    # The variance of a sum of 2k independent Rademacher variables is 2k.
    var_X = 2 * k

    # The components of the sum S are:
    # Sx = 1*X_A + 0.5*X_B - 0.5*X_C
    # Sy = (sqrt(3)/2)*X_B + (sqrt(3)/2)*X_C
    
    # The variances of Sx and Sy are calculated based on the independence of X_A, X_B, X_C.
    var_Sx = 1**2 * var_X + (sp.S(1)/2)**2 * var_X + (-sp.S(1)/2)**2 * var_X
    var_Sy = (sp.sqrt(3)/2)**2 * var_X + (sp.sqrt(3)/2)**2 * var_X
    
    # The covariance Cov(Sx, Sy) is 0.
    # Cov(Sx, Sy) = Cov(X_A + 0.5*X_B - 0.5*X_C, (sqrt(3)/2)*(X_B + X_C))
    #            = (sqrt(3)/2) * [0.5*Var(X_B) - 0.5*Var(X_C)] = 0.

    # Simplify the variance expressions. They are both equal.
    var_S = sp.simplify(var_Sx)
    
    print(f"\nThe variance of both Sx and Sy is found to be 3*k. So, sigma^2 = {var_S}")

    # Step 2: Express the variance in terms of n.
    # Since n = 6k, we have k = n/6.
    var_S_n = var_S.subs(k, n/6)
    print(f"Substituting k = n/6, the variance becomes sigma^2 = {var_S_n}")

    # Step 3: Calculate the probability P(n) = P(||S||_2^2 <= 2).
    # For a 2D Gaussian with independent components of variance sigma^2,
    # the probability P(Sx^2 + Sy^2 <= R^2) is given by 1 - exp(-R^2 / (2 * sigma^2)).
    R_squared = 2
    P_n_expr = 1 - sp.exp(-R_squared / (2 * var_S_n))
    P_n_simplified = sp.simplify(P_n_expr)
    
    print("\nStep 2: Calculating the probability P(n) = P(||S||_2 <= sqrt(2)).")
    print(f"The formula for P(n) for large n is: P(n) = {P_n_simplified}")

    # Step 4: Calculate the limit of n * P(n) as n -> infinity.
    expression_to_limit = n * P_n_simplified
    
    print("\nStep 3: Computing the final limit of n*P(n) as n approaches infinity.")
    
    # Create a string representation of the final equation
    # The numbers in the equation are n, 1, -2, n
    final_equation_str = f"lim (n * (1 - exp(-2/n))) as n -> oo"
    print(f"We need to compute: {final_equation_str}")

    limit_value = sp.limit(expression_to_limit, n, sp.oo)

    print(f"\nThe result of the limit is: {limit_value}")
    
    return limit_value

if __name__ == '__main__':
    solve_limit_problem()