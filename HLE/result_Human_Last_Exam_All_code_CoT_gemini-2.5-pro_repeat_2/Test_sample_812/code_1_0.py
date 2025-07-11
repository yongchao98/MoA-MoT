import sympy
from sympy import Abs, oo, integrate, Symbol, Rational, diff

def solve_cumulant():
    """
    This function calculates the third cumulant of the converged variable Y_n.
    """
    # Define symbolic variables
    x = Symbol('x')
    t = Symbol('t')

    # Define the PDF of the random variables X_i
    pdf = Rational(3, 2) / (1 + Abs(x))**4

    print("This script calculates the third cumulant of a limiting distribution based on the Central Limit Theorem.")
    print("It assumes the standard formulation Y_n = sqrt(n) * (X_bar - mu), where X_bar is the sample mean.\n")
    
    print("--- Step 1: Analyze the distribution of X_i ---")
    print(f"The probability density function is f(x) = {pdf}\n")

    # --- Calculate Mean (mu) ---
    print("--- Step 2: Calculate the mean (mu) of X_i ---")
    # Due to the symmetry of f(x) around 0, the integral of x*f(x) is 0.
    mu = integrate(x * pdf, (x, -oo, oo))
    print(f"The mean is mu = E[X] = integral from -inf to inf of x*f(x) dx.")
    print(f"Calculated mean: mu = {mu}\n")

    # --- Calculate Variance (sigma^2) ---
    print("--- Step 3: Calculate the variance (sigma^2) of X_i ---")
    # E[X^2]
    E_x2 = integrate(x**2 * pdf, (x, -oo, oo))
    # Variance sigma^2 = E[X^2] - mu^2
    sigma2 = E_x2 - mu**2
    print(f"The variance is sigma^2 = E[X^2] - mu^2.")
    print(f"First, we calculate E[X^2] = integral from -inf to inf of x^2*f(x) dx = {E_x2}")
    print(f"Then, the variance is sigma^2 = {E_x2} - ({mu})^2 = {sigma2}\n")

    # --- Apply Central Limit Theorem ---
    print("--- Step 4: Apply the Central Limit Theorem ---")
    print("As n -> infinity, Y_n converges in distribution to a Normal random variable Y.")
    print(f"The distribution of Y is N(mean=0, variance=sigma^2), which is N(0, {sigma2}).\n")

    # --- Calculate Third Cumulant ---
    print("--- Step 5: Calculate the third cumulant of the limiting Normal distribution Y ---")
    # The cumulant-generating function (CGF) for a Normal(mean, variance) is K(t) = mean*t + variance*t^2/2
    limit_mean = 0
    limit_variance = sigma2
    K = limit_mean * t + limit_variance * t**2 / 2
    
    print(f"The CGF of Y ~ N({limit_mean}, {limit_variance}) is K(t) = {limit_mean}*t + {limit_variance}*t^2/2.")
    print(f"The equation for K(t) is: K(t) = {K}")

    # The k-th cumulant is the k-th derivative of K(t) with respect to t, evaluated at t=0.
    K_d1 = diff(K, t)
    K_d2 = diff(K_d1, t)
    K_d3 = diff(K_d2, t)

    third_cumulant = K_d3.subs(t, 0)
    
    print("To find the third cumulant, we calculate the third derivative of K(t) and evaluate it at t=0:")
    print(f"dK/dt = {K_d1}")
    print(f"d^2K/dt^2 = {K_d2}")
    print(f"d^3K/dt^3 = {K_d3}")
    print(f"The third cumulant kappa_3 = (d^3K/dt^3 at t=0) = {third_cumulant}\n")

    print("--- Final Answer ---")
    print(f"The third cumulant of the converged variable is: {third_cumulant}")

if __name__ == "__main__":
    solve_cumulant()
    print("\n<<<0>>>")
