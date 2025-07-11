import sympy as sp

def solve_probability_limit():
    """
    This function calculates the limit of n*P(n) as n -> infinity using symbolic mathematics.
    """
    # Define symbols for our calculation
    # n is the total number of vectors, which is 6k.
    n = sp.Symbol('n', positive=True)

    # --- Plan Step 1 & 2: Analyze the sum S and apply CLT ---
    # The sum is S = (X, Y). We need the covariance matrix of (X,Y).
    # We have 3 groups of vectors: 2k of v_A, 2k of v_B, 2k of v_C.
    # Total vectors n = 6k, so each group has n/3 vectors.
    # Let S_A, S_B, S_C be the sum of Rademacher variables for each group.
    # E.g., S_A = sum_{i=1}^{n/3} epsilon_i.
    # By CLT, S_A, S_B, S_C are approximately N(0, n/3).
    var_S_group = n / 3

    # The components of the sum S are:
    # X = 1*S_A + 0.5*S_B - 0.5*S_C
    # Y = 0*S_A + (sqrt(3)/2)*S_B + (sqrt(3)/2)*S_C
    # Since S_A, S_B, S_C are independent, we can calculate Var(X) and Var(Y).

    print("Step 1: Calculating the covariance matrix of the sum vector S=(X,Y)")
    
    # --- Plan Step 3: Calculate the covariance matrix ---
    var_X = (1**2) * var_S_group + (sp.Rational(1, 2))**2 * var_S_group + (-sp.Rational(1, 2))**2 * var_S_group
    var_Y = (sp.sqrt(3)/2)**2 * var_S_group + (sp.sqrt(3)/2)**2 * var_S_group
    
    # The covariance term uses independence of S_A, S_B, S_C.
    # Cov(X,Y) = Cov(S_A + 0.5*S_B - 0.5*S_C, (sqrt(3)/2)*S_B + (sqrt(3)/2)*S_C)
    # Due to independence, only terms with the same S-variable are non-zero.
    # Cov(X,Y) = 0.5*(sqrt(3)/2)*Var(S_B) - 0.5*(sqrt(3)/2)*Var(S_C)
    cov_XY = sp.Rational(1, 2) * (sp.sqrt(3)/2) * var_S_group - sp.Rational(1, 2) * (sp.sqrt(3)/2) * var_S_group
    
    var_X = sp.simplify(var_X)
    var_Y = sp.simplify(var_Y)

    print(f"Variance of X: Var(X) = {var_X}")
    print(f"Variance of Y: Var(Y) = {var_Y}")
    print(f"Covariance of X and Y: Cov(X, Y) = {cov_XY}")
    print("The covariance matrix is diagonal, Sigma = [[n/2, 0], [0, n/2]].")
    print("-" * 30)

    # --- Plan Step 4: Calculate the probability P(n) ---
    # The PDF of S=(X,Y) is a bivariate normal distribution f(x,y).
    # f(x,y) = 1/(2*pi*det(Sigma)^0.5) * exp(-0.5 * (x,y) * Sigma_inv * (x,y)^T)
    # With Sigma = [[n/2, 0], [0, n/2]], this simplifies to:
    # f(x,y) = 1/(pi*n) * exp(-(x^2+y^2)/n)
    # We need to integrate this over the disk x^2+y^2 <= 2.
    
    print("Step 2: Calculating the probability P(n)")
    print("P(n) is the integral of the PDF over the disk with radius sqrt(2).")

    r, theta = sp.symbols('r theta', positive=True)
    # PDF in polar coordinates
    pdf_polar = 1 / (sp.pi * n) * sp.exp(-r**2 / n)
    # Integrand with Jacobian r
    integrand = pdf_polar * r
    
    # Integrate over r from 0 to sqrt(2)
    integral_r = sp.integrate(integrand, (r, 0, sp.sqrt(2)))
    # Integrate over theta from 0 to 2*pi
    P_n = sp.integrate(integral_r, (theta, 0, 2 * sp.pi))

    print(f"The probability P(n) is approximately: P(n) = {P_n}")
    print("-" * 30)
    
    # --- Plan Step 5: Evaluate the Limit ---
    print("Step 3: Calculating the final limit of n*P(n)")
    
    # The expression to take the limit of
    limit_expr = n * P_n
    
    # Calculate the limit as n -> infinity
    result = sp.limit(limit_expr, n, sp.oo)
    
    print(f"The expression is: n * ({P_n})")
    print(f"The limit as n approaches infinity is: {result}")

if __name__ == '__main__':
    solve_probability_limit()