import math

def estimate_fluctuation_magnitude(epsilon):
    """
    Estimates the maximum magnitude of fluctuations R(epsilon) based on the derived formula.

    The formula is R(epsilon) = sqrt(epsilon / 192).

    Args:
        epsilon (float): A small positive parameter, 0 < epsilon < 1.

    Returns:
        float: The estimated value of R.
    """
    if not (0 < epsilon < 1):
        print("Error: Epsilon must be between 0 and 1.")
        return None

    # The formula for R^2 is epsilon / 192
    r_squared_numerator = epsilon
    r_squared_denominator = 192
    r_squared = r_squared_numerator / r_squared_denominator

    # R is the square root of R^2
    R = math.sqrt(r_squared)

    print(f"For epsilon = {epsilon}:")
    print(f"  The scaling relationship is R(epsilon) ~ epsilon^(1/2)")
    print(f"  The equation for R^2 is: R^2 = epsilon / 192")
    print(f"  Calculating the values:")
    print(f"    Numerator (epsilon): {r_squared_numerator}")
    print(f"    Denominator: {r_squared_denominator}")
    print(f"  R^2 = {r_squared_numerator} / {r_squared_denominator} = {r_squared:.6e}")
    print(f"  R = sqrt(R^2) = {R:.6e}")
    print("-" * 20)
    
    return R

# --- Main execution ---
# You can test with different values of epsilon
epsilon_values = [0.1, 0.05, 0.01]
for eps in epsilon_values:
    estimate_fluctuation_magnitude(eps)
