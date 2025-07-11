import math

def calculate_it3_mf_bound(x, u, c, sigma):
    """
    Calculates the membership grade for a Gaussian-based Interval Type-3 MF bound.

    This formulation represents a bound (e.g., lower bound) of the IT3 MF,
    where the uncertainty is parameterized by the secondary variable 'u', which
    scales the height of a Gaussian function.

    Args:
        x (float): The primary input variable.
        u (float): The secondary input variable (in [0, 1]) for uncertainty.
        c (float): The center (mean) of the Gaussian function.
        sigma (float): The standard deviation (width) of the Gaussian function.

    Returns:
        float: The calculated membership grade for the specified bound.
    """
    if not 0 <= u <= 1:
        raise ValueError("Secondary variable 'u' must be in the interval [0, 1].")
    if sigma <= 0:
        raise ValueError("Standard deviation 'sigma' must be positive.")

    # Calculate the exponent part of the Gaussian function
    exponent = -0.5 * ((x - c) / sigma)**2
    
    # Calculate the value of the Gaussian function
    gaussian_value = math.exp(exponent)
    
    # Scale the Gaussian by the secondary variable 'u'
    membership_grade = u * gaussian_value
    
    return membership_grade, exponent, gaussian_value

# --- Example Parameters ---
# The primary input value
x_input = 6.0
# The secondary variable, representing the uncertainty level (between 0 and 1)
u_input = 0.9
# The center of the Gaussian membership function
center = 5.0
# The standard deviation (width) of the Gaussian membership function
std_dev = 2.0

# --- Calculation ---
try:
    # Calculate the membership grade
    mu_result, exp_val, gauss_val = calculate_it3_mf_bound(x_input, u_input, center, std_dev)

    # --- Output ---
    # Print the general mathematical formulation
    print("The mathematical formulation for the IT3 MF bound is:")
    print("μ_bound(x, u) = u * exp(-0.5 * ((x - c) / σ)^2)\n")

    # Print the equation with the specific values used in the calculation
    print("Substituting the example values:")
    print(f"x = {x_input}, u = {u_input}, c = {center}, σ = {std_dev}\n")
    
    print("The final equation with calculated values is:")
    final_equation = f"μ_bound({x_input}, {u_input}) = {u_input} * exp(-0.5 * (({x_input} - {center}) / {std_dev})^2) = {mu_result:.6f}"
    print(final_equation)

except ValueError as e:
    print(f"Error: {e}")
