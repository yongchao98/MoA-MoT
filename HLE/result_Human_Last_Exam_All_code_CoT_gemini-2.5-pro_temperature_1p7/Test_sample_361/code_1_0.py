import math

def calculate_critical_potential(lambda_gb, L=1.0):
    """
    Calculates the critical chemical potential for scalar condensation in a holographic model
    with Einstein-Gauss-Bonnet gravity, based on an analytic formula from the D4/D8 model.

    Args:
        lambda_gb (float): The Gauss-Bonnet coupling constant.
        L (float): The AdS radius.

    Returns:
        None. Prints the calculation and the final result.
    """

    # Calculate the term inside the square root of L_eff_sq definition
    # This term must be non-negative, which constrains lambda_gb <= 1/4
    if 1.0 - 4.0 * lambda_gb < 0:
        print("Error: The Gauss-Bonnet coupling lambda_gb must be <= 1/4.")
        return

    # Calculate the effective AdS radius squared, L_eff_sq
    one_minus_sqrt_term = 1.0 - math.sqrt(1.0 - 4.0 * lambda_gb)
    L_eff_sq_numerator = L**2 * one_minus_sqrt_term
    L_eff_sq_denominator = 2.0 * lambda_gb
    L_eff_sq = L_eff_sq_numerator / L_eff_sq_denominator

    # Calculate the term involving the Gamma function
    gamma_quarter_pow4 = math.gamma(0.25)**4

    # Calculate the constant denominator
    denominator = 32.0 * math.pi**2

    # Calculate the square of the critical chemical potential
    mu_c_squared = L_eff_sq * gamma_quarter_pow4 / denominator

    # Calculate the final critical chemical potential
    mu_c = math.sqrt(mu_c_squared)

    # Print the full calculation, showing each number in the final equation
    print(f"sqrt(( ({L**2:.1f} * (1 - sqrt(1 - 4 * {lambda_gb}))) / (2 * {lambda_gb}) ) * {gamma_quarter_pow4} / {denominator}) = {mu_c}")


# Set the given Gauss-Bonnet coupling
gauss_bonnet_coupling = 0.1

# Run the calculation and print the result
calculate_critical_potential(gauss_bonnet_coupling)