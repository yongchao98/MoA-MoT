import numpy as np
from scipy.integrate import quad

def solve():
    """
    This script calculates the third cumulant of the limiting distribution
    of the random variable Y_n specified in the problem.
    """

    # 1. Define the probability density function f(x)
    # f(x) = 3 / (2 * (1 + |x|)^4)
    def f(x):
        return (3/2) / (1 + np.abs(x))**4

    # 2. Verify the conditions for the Central Limit Theorem (finite mean and variance)
    # by calculating the first few moments of X.

    # Integrand for the mean E[X] = integral(x * f(x) dx)
    def x_f(x):
        return x * f(x)

    # Integrand for the second raw moment E[X^2] = integral(x^2 * f(x) dx)
    def x2_f(x):
        return x**2 * f(x)
        
    print("Step 1: Calculate the moments of the X_i variables to verify CLT conditions.")
    
    # Due to the symmetry of f(x) around 0, odd moments like the mean are zero.
    # The integral of an odd function over a symmetric interval is 0.
    # E[X] = integral from -inf to inf of x * f(x) dx = 0
    mu, mu_err = quad(x_f, -np.inf, np.inf)
    print(f"The mean (mu) of X_i is: {mu:.6f}")

    # The variance Var(X) = E[X^2] - (E[X])^2. Since mu=0, Var(X) = E[X^2].
    mu2_prime, mu2_err = quad(x2_f, -np.inf, np.inf)
    variance = mu2_prime - mu**2
    print(f"The variance (sigma^2) of X_i is: {variance:.6f}")
    
    # 3. Apply the Central Limit Theorem
    print("\nStep 2: Apply the Central Limit Theorem (CLT).")
    print("Since the mean and variance are finite, the CLT applies.")
    print("The variable Y_n = sqrt(n) * (sample_mean(X) - mu) converges in distribution to a Normal distribution.")
    print(f"The limiting distribution is N(0, sigma^2), which is N(0, {variance:.1f}).")
    
    # 4. Determine the third cumulant of the limiting Normal distribution.
    print("\nStep 3: Find the third cumulant of the resulting Normal distribution.")
    # The cumulants of a Normal distribution N(mu_norm, sigma_norm^2) are:
    # k_1 = mu_norm, k_2 = sigma_norm^2, and k_j = 0 for j >= 3.
    # For our limiting distribution N(0, 1.0), the parameters are mu_norm=0 and sigma_norm^2=1.0.
    third_cumulant_of_limit = 0
    
    print("The third cumulant of any Normal distribution is always 0.")
    print(f"The equation for the third cumulant of the converged variable is simply:")
    print(f"{third_cumulant_of_limit}")

solve()