import math

def estimate_fluctuation_magnitude(epsilon):
    """
    Estimates the maximum magnitude of fluctuations (R) as a function of epsilon.

    Args:
        epsilon (float): A small parameter, 0 < epsilon < 1.

    Returns:
        float: The estimated value of R.
    """
    if not (0 < epsilon < 1):
        print("Error: epsilon must be between 0 and 1.")
        return None

    # From the analytical derivation, the maximum variance is given by:
    # max_Var = epsilon / 192
    # R is the square root of the maximum variance.
    # R^2 = epsilon / 192
    
    constant_C_squared = 1 / 192
    max_variance = epsilon * constant_C_squared
    R = math.sqrt(max_variance)
    
    print(f"For epsilon = {epsilon}:")
    print("The maximum variance of the fluctuations is estimated as:")
    print(f"Var_max = epsilon / 192 = {epsilon} / 192 = {max_variance}")
    print("\nThe magnitude R is the square root of the maximum variance:")
    print(f"R = (epsilon / 192)^(1/2) = ({epsilon} / 192)^(1/2)")
    print(f"The calculated value is:")
    print(R)
    
    return R

# --- Main execution ---
# Set a value for epsilon
epsilon_value = 0.01

# Calculate and print the result
estimate_fluctuation_magnitude(epsilon_value)

<<<0.07216878364870322 * epsilon**(0.5)>>>