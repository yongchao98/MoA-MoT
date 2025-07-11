import numpy as np
from scipy.integrate import quad

def solve_cumulant():
    """
    This function calculates the third cumulant of the converged variable Yn.
    """
    
    # Define the probability density function (PDF) of X
    # f(x) = 3 / (2 * (1 + |x|)^4)
    
    # Step 1: Calculate the mean (mu) of X.
    # mu = E[X] = integral from -inf to inf of x * f(x) dx
    # The integrand x * f(x) is an odd function, so the integral over a symmetric domain is 0.
    mean_integrand = lambda x: x * 3 / (2 * (1 + abs(x))**4)
    mu, _ = quad(mean_integrand, -np.inf, np.inf)
    
    # Step 2: Calculate the second moment E[X^2] of X.
    # E[X^2] = integral from -inf to inf of x^2 * f(x) dx
    second_moment_integrand = lambda x: x**2 * 3 / (2 * (1 + abs(x))**4)
    e_x2, _ = quad(second_moment_integrand, -np.inf, np.inf)
    
    # Step 3: Calculate the variance (sigma^2) of X.
    # sigma^2 = E[X^2] - mu^2
    variance = e_x2 - mu**2
    
    print("Step 1: Calculate the mean and variance of the random variables X_i.")
    print(f"The mean E[X] is: {mu:.4f}")
    print(f"The second moment E[X^2] is: {e_x2:.4f}")
    print(f"The variance Var(X) = E[X^2] - (E[X])^2 = {e_x2:.4f} - {mu:.4f}^2 = {variance:.4f}")
    print("-" * 30)

    # Step 4: Apply the Central Limit Theorem.
    # The variable Y_n converges in distribution to a Normal distribution N(0, sigma^2).
    print("Step 2: Identify the limiting distribution using the Central Limit Theorem.")
    print(f"Y_n converges to a Normal distribution with mean 0 and variance {variance:.4f}, i.e., N(0, {variance:.4f}).")
    print("-" * 30)
    
    # Step 5: Find the third cumulant of the limiting distribution.
    # For any Normal distribution, the cumulants of order k >= 3 are all zero.
    third_cumulant = 0
    
    print("Step 3: Find the third cumulant of the limiting Normal distribution.")
    print("For any Normal distribution, the third cumulant is 0.")
    print(f"The final equation for the third cumulant (k3) is: k3 = {third_cumulant}")
    print("-" * 30)

    print(f"The third cumulant of the converged variable is: {third_cumulant}")

solve_cumulant()