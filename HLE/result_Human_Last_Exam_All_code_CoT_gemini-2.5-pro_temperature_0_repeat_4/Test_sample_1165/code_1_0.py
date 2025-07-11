import math

def estimate_fluctuations(epsilon):
    """
    Estimates the maximum magnitude of fluctuations R as a function of epsilon.

    Args:
        epsilon (float): A small parameter, 0 < epsilon < 1.
    """
    if not (0 < epsilon < 1):
        print("Error: epsilon must be between 0 and 1.")
        return

    # The derived formula for the maximum variance of the second-order term y2(x)
    # max_x Var[y2(x)] = 1 / (192 * epsilon^3)
    c = 192
    max_var_y2_numerator = 1
    max_var_y2_denominator = c * (epsilon**3)
    max_var_y2 = max_var_y2_numerator / max_var_y2_denominator

    # The formula for R^2 is epsilon^4 * max_x Var[y2(x)]
    r_squared_factor = epsilon**4
    r_squared = r_squared_factor * max_var_y2

    # R is the square root of R^2
    r = math.sqrt(r_squared)

    print(f"For epsilon = {epsilon}:")
    print("-" * 30)
    print("Step 1: Calculate the maximum variance of y2(x)")
    print(f"max_x Var[y2(x)] = {max_var_y2_numerator} / ({c} * {epsilon}^3)")
    print(f"                   = {max_var_y2_numerator} / ({c} * {epsilon**3:.4g})")
    print(f"                   = {max_var_y2_numerator} / {max_var_y2_denominator:.4g}")
    print(f"                   = {max_var_y2:.4g}")
    print()
    print("Step 2: Calculate R^2")
    print(f"R^2 = epsilon^4 * max_x Var[y2(x)]")
    print(f"    = {epsilon}^4 * {max_var_y2:.4g}")
    print(f"    = {r_squared_factor:.4g} * {max_var_y2:.4g}")
    print(f"    = {r_squared:.4g}")
    print()
    print("Step 3: Calculate R")
    print(f"R = sqrt(R^2) = sqrt({r_squared:.4g})")
    print(f"  = {r:.4g}")
    print("-" * 30)
    print("The final estimated value for R is:", r)

# Example calculation for epsilon = 0.1
estimate_fluctuations(0.1)
