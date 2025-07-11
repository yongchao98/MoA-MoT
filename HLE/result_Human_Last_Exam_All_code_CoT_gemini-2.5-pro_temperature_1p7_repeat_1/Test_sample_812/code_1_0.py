import sympy

def solve_cumulant_problem():
    """
    Calculates the mean and variance of the given distribution to verify
    the conditions for the Central Limit Theorem, and then determines the
    third cumulant of the resulting limiting Normal distribution.
    """

    # Step 1: Define the variables and the probability density function (PDF)
    x = sympy.Symbol('x', real=True)
    # The PDF is f(x) = 3 / (2 * (1 + abs(x))^4).
    # To handle the absolute value, we split the function for x > 0 and x < 0.
    pdf_pos = 3 / (2 * (1 + x)**4)  # for x >= 0
    pdf_neg = 3 / (2 * (1 - x)**4)  # for x < 0

    print("This script solves for the third cumulant of the converged variable.")
    print("Plan:")
    print("1. Verify conditions for the Central Limit Theorem (CLT) by calculating the mean and variance of X_i.")
    print("2. Identify the limiting distribution based on the CLT.")
    print("3. State the third cumulant of this limiting distribution.")
    print("-" * 50)

    # Step 2: Verify CLT conditions by calculating mean and variance
    print("Step 1: Calculating Mean (mu) and Variance (sigma^2)")
    # Calculate the mean (mu = E[X])
    # E[X] = integral from -inf to inf of x*f(x) dx
    print("\nCalculating the mean of X_i, mu = E[X]:")
    integral_mean_neg = sympy.integrate(x * pdf_neg, (x, -sympy.oo, 0))
    integral_mean_pos = sympy.integrate(x * pdf_pos, (x, 0, sympy.oo))
    mu = integral_mean_neg + integral_mean_pos

    print(f"  Integral over (-inf, 0): {integral_mean_neg}")
    print(f"  Integral over (0, inf):  {integral_mean_pos}")
    print(f"  Total Mean mu = {integral_mean_neg} + {integral_mean_pos} = {mu}")
    print("  Conclusion: The mean is finite.")

    # Calculate the variance (sigma^2 = E[X^2] - mu^2)
    # First, calculate the second moment E[X^2] = integral of x^2*f(x) dx
    print("\nCalculating the variance of X_i, sigma^2 = E[X^2] - mu^2:")
    integral_E_x2_neg = sympy.integrate(x**2 * pdf_neg, (x, -sympy.oo, 0))
    integral_E_x2_pos = sympy.integrate(x**2 * pdf_pos, (x, 0, sympy.oo))
    E_x2 = integral_E_x2_neg + integral_E_x2_pos
    variance = E_x2 - mu**2
    
    print(f"  Calculating the second moment, E[X^2]:")
    print(f"    Integral over (-inf, 0): {integral_E_x2_neg}")
    print(f"    Integral over (0, inf):  {integral_E_x2_pos}")
    print(f"    Total Second Moment E[X^2] = {integral_E_x2_neg} + {integral_E_x2_pos} = {E_x2}")
    print(f"  Variance sigma^2 = E[X^2] - mu^2 = {E_x2} - ({mu})**2 = {variance}")
    print("  Conclusion: The variance is finite and non-zero.")
    print("-" * 50)

    # Step 3: Apply CLT and find the third cumulant
    print("Step 2: Identifying the Limiting Distribution")
    print("Since the mean and variance are finite, the Central Limit Theorem applies.")
    print(f"The variable Y_n converges in distribution to a Normal random variable Y ~ N(0, sigma^2).")
    print(f"In this case, the limiting distribution is Y ~ Normal(mean=0, variance={variance}).")
    print("-" * 50)
    
    print("Step 3: Finding the Third Cumulant")
    third_cumulant_normal = 0
    print("A fundamental property of any Normal distribution is that all its cumulants of order 3 or higher are zero.")
    print(f"The first cumulant is the mean (0).")
    print(f"The second cumulant is the variance ({variance}).")
    print(f"The third cumulant is {third_cumulant_normal}.")

solve_cumulant_problem()