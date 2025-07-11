import sympy

def solve_limit():
    """
    This function calculates the limit of n*P(n) as n -> infinity.
    The derivation is based on the Central Limit Theorem.
    """
    # Define k as a symbolic variable representing a large positive integer.
    # As k -> infinity, n -> infinity.
    k = sympy.Symbol('k', positive=True)

    # The total number of vectors is n = 6k.
    n = 6 * k
    print(f"The total number of vectors, n, is expressed in terms of k as: {n}")

    # The sum S = (Sx, Sy) is approximately a 2D Gaussian vector N(0, Sigma).
    # From the problem statement, we calculate the variance of each component of S.
    # Var(Sx) = Var(Sy) = 3k.
    # The components are uncorrelated, so the covariance matrix is diagonal.
    variance_per_component = 3 * k

    # We are interested in the probability that ||S|| <= sqrt(2), which means ||S||^2 <= 2.
    R_squared = 2

    # For a 2D isotropic Gaussian distribution with variance sigma^2 in each component,
    # the probability P(||S||^2 <= R^2) is given by 1 - exp(-R^2 / (2 * sigma^2)).
    # Here, sigma^2 corresponds to our variance_per_component.
    prob_n = 1 - sympy.exp(-R_squared / (2 * variance_per_component))
    print(f"The probability P(n) for large n is approximated by the expression: {prob_n}")

    # We need to find the limit of the product n * P(n).
    expression_to_limit = n * prob_n
    print(f"The expression for n*P(n) to be evaluated in the limit is: {expression_to_limit}")

    # Calculate the limit as k -> oo.
    limit_result = sympy.limit(expression_to_limit, k, sympy.oo)

    print(f"\nThe limit of n*P(n) as n goes to infinity is:")
    print(limit_result)

solve_limit()