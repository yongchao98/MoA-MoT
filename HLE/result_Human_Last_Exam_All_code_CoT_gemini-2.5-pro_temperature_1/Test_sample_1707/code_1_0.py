import decimal

def solve_for_T():
    """
    Calculates the time T based on the solvability condition of the boundary-value problem.
    """
    # The numbers involved have very large and small exponents, so we need high precision.
    # The exponent of the result is around 5,010,000. We set the precision
    # of the decimal context to a slightly higher value to ensure accuracy.
    decimal.getcontext().prec = 5010010

    # Define the given parameters as Decimal objects for high-precision arithmetic.
    alpha = decimal.Decimal('1e10000')
    x0 = decimal.Decimal('1e-5000000')
    y0 = decimal.Decimal('1e-5000000')

    # From the perturbation analysis, the solvability condition is KT = alpha,
    # where K = x0/(1-x0) + y0/(1-y0).
    # We solve for T = alpha / K.

    # Calculate K
    K = x0 / (1 - x0) + y0 / (1 - y0)

    # Calculate T
    T = alpha / K

    # The problem requires outputting each number in the final equation, which is T = alpha / K.
    # We will print the values of alpha, K, and the resulting T.
    # We use scientific notation for readability.
    
    print(f"The final equation to be solved is T = alpha / K")
    print(f"where K = x0/(1-x0) + y0/(1-y0).")
    print(f"Using the given values:")
    print(f"alpha = {alpha:.1e}")
    print(f"x0 = {x0:.1e}")
    print(f"y0 = {y0:.1e}")
    print("-" * 20)
    print(f"Calculated intermediate value K = {K:.6e}")
    print(f"Final calculation: T = {alpha:.1e} / {K:.6e}")
    print(f"Resulting time T = {T:.6e}")

solve_for_T()