import math

def calculate_it3_gaussian_bounds(x, u, c, sigma_l, sigma_u):
    """
    Calculates the vertical cross-section bounds of a Gaussian IT3 MF.

    This function implements a common model where the standard deviation
    of the Gaussian is modulated by the secondary variable 'u'.

    Args:
        x (float): The primary input variable.
        u (float): The secondary input variable, in [0, 1].
        c (float): The center of the Gaussian function.
        sigma_l (float): The lower bound of the primary standard deviation range.
        sigma_u (float): The upper bound of the primary standard deviation range.

    Returns:
        tuple: A tuple containing (mu_lower, mu_upper).
    """
    if not (0 <= u <= 1):
        raise ValueError("Secondary variable u must be in the range [0, 1].")
    if sigma_l > sigma_u:
        raise ValueError("sigma_l cannot be greater than sigma_u.")
    if sigma_l <= 0 or sigma_u <= 0:
        raise ValueError("Standard deviations must be positive.")

    # Calculate the average and radius of the standard deviation interval
    sigma_avg = (sigma_l + sigma_u) / 2.0
    sigma_rad = (sigma_u - sigma_l) / 2.0

    # Calculate the u-dependent standard deviations
    # sigma_lower_u corresponds to the Lower Membership Function (LMF)
    # sigma_upper_u corresponds to the Upper Membership Function (UMF)
    sigma_lower_u = sigma_avg - u * sigma_rad
    sigma_upper_u = sigma_avg + u * sigma_rad

    # Calculate the exponent term (x - c)^2
    x_term_sq = (x - c)**2

    # The Gaussian is an increasing function of sigma.
    # Therefore, the smaller sigma_lower_u yields the lower membership bound.
    exponent_lower = -0.5 * (x_term_sq / (sigma_lower_u**2))
    mu_lower = math.exp(exponent_lower)

    # The larger sigma_upper_u yields the upper membership bound.
    exponent_upper = -0.5 * (x_term_sq / (sigma_upper_u**2))
    mu_upper = math.exp(exponent_upper)

    return mu_lower, mu_upper, sigma_lower_u, sigma_upper_u


# --- Example Calculation ---

# Define the parameters for the IT3 MF
x_prime = 5.5  # A fixed value for the primary variable x
u_val = 0.8    # A fixed value for the secondary variable u
center = 5.0   # Center of the Gaussian
sigma_lower_bound = 1.0 # Lower bound of the std. dev. range
sigma_upper_bound = 2.0 # Upper bound of the std. dev. range

# Calculate the bounds
try:
    (
        mu_l,
        mu_u,
        s_l_u,
        s_u_u,
    ) = calculate_it3_gaussian_bounds(
        x_prime, u_val, center, sigma_lower_bound, sigma_upper_bound
    )

    # --- Print the formulation and the results ---
    print("Formulation for the Vertical Cross-Section of a Gaussian IT3 MF:")
    print("-" * 65)
    print("μ_lower(x,u) = exp( -0.5 * ( (x - c) / σ_L(u) )^2 )")
    print("μ_upper(x,u) = exp( -0.5 * ( (x - c) / σ_U(u) )^2 )")
    print("\nwhere:")
    print("σ_L(u) = (σ_lower + σ_upper)/2 - u * (σ_upper - σ_lower)/2")
    print("σ_U(u) = (σ_lower + σ_upper)/2 + u * (σ_upper - σ_lower)/2")
    print("-" * 65)

    print("\nCalculation with the given values:")
    print(f"x = {x_prime}, u = {u_val}, c = {center}, σ_lower = {sigma_lower_bound}, σ_upper = {sigma_upper_bound}\n")
    
    # Print the derivation for sigma_L(u) and sigma_U(u)
    sigma_avg = (sigma_lower_bound + sigma_upper_bound) / 2.0
    sigma_rad = (sigma_upper_bound - sigma_lower_bound) / 2.0
    print(f"σ_L({u_val}) = ({sigma_lower_bound} + {sigma_upper_bound})/2 - {u_val} * ({sigma_upper_bound} - {sigma_lower_bound})/2 = {sigma_avg} - {u_val} * {sigma_rad} = {s_l_u:.4f}")
    print(f"σ_U({u_val}) = ({sigma_lower_bound} + {sigma_upper_bound})/2 + {u_val} * ({sigma_upper_bound} - {sigma_lower_bound})/2 = {sigma_avg} + {u_val} * {sigma_rad} = {s_u_u:.4f}\n")
    
    # Print the final equations with all numbers
    print("Resulting Bounds:")
    print(f"μ_lower = exp( -0.5 * ( ({x_prime} - {center}) / {s_l_u:.4f} )^2 ) = {mu_l:.4f}")
    print(f"μ_upper = exp( -0.5 * ( ({x_prime} - {center}) / {s_u_u:.4f} )^2 ) = {mu_u:.4f}\n")
    
    print(f"The vertical cross-section at x={x_prime} and u={u_val} is the interval: [{mu_l:.4f}, {mu_u:.4f}]")

except ValueError as e:
    print(f"Error: {e}")
