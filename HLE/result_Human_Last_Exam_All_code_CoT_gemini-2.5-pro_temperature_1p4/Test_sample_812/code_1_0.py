import sympy
from sympy import symbols, Abs, integrate, oo, Rational

def solve_cumulant_problem():
    """
    This function solves the problem of finding the third cumulant of the
    converged random variable Y_n by following these steps:
    1. Clarify the definition of Y_n based on the Central Limit Theorem (CLT).
    2. Calculate the mean (μ) and variance (σ^2) of X_i from its PDF to check CLT conditions.
    3. Identify the limit distribution of Y_n.
    4. Find the third cumulant of this limit distribution.
    """
    # Define the symbolic variable and the probability density function (PDF)
    x = symbols('x')
    pdf = Rational(3, 2) / (1 + Abs(x))**4

    print("Step 1: Calculating the mean (μ) and variance (σ^2) of X_i.")
    
    # Calculate the mean (μ) by integrating x * pdf from -∞ to +∞.
    # Sympy's integrate function can handle the absolute value.
    mu = integrate(x * pdf, (x, -oo, oo))
    
    # Calculate the second moment E[X^2] by integrating x^2 * pdf.
    e_x_squared = integrate(x**2 * pdf, (x, -oo, oo))
    
    # Calculate the variance (σ^2) using the formula σ^2 = E[X^2] - μ^2.
    sigma_squared = e_x_squared - mu**2

    print(f"   - The mean of X_i is μ = {mu}.")
    print(f"   - The variance of X_i is σ^2 = {sigma_squared}.")

    print("\nStep 2: Applying the Central Limit Theorem (CLT).")
    print("   The random variable Y_n is assumed to be the standardized sample mean,")
    print("   Y_n = sqrt(n) * (mean(X_i) - μ), which converges in distribution.")
    print(f"   As n → ∞, Y_n converges to a Normal distribution Y ~ N(0, σ^2).")
    print(f"   Given our calculated variance, the limit distribution is Y ~ N(0, {sigma_squared}).")

    print("\nStep 3: Determining the third cumulant of the limit distribution.")
    print("   For any Normal distribution, all cumulants of order 3 and higher are 0.")
    
    # The third cumulant of a Normal distribution is always 0.
    third_cumulant_of_limit = 0
    
    print("\nFinal Result:")
    print("The final equation for the third cumulant (κ₃) of the converged variable is:")
    print(f"κ₃ = {third_cumulant_of_limit}")


# Execute the solver function
solve_cumulant_problem()