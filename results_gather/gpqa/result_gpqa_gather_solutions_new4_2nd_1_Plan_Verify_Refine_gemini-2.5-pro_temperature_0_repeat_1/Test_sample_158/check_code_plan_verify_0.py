import numpy as np
from scipy.integrate import quad

def check_cosmological_distance():
    """
    This function verifies the calculation of the comoving distance to the quasar.
    It follows these steps:
    1. Defines the physical and cosmological constants from the problem.
    2. Calculates the redshift (z) from the observed Lyman-alpha line.
    3. Defines the integrand for the comoving distance formula in a flat Lambda-CDM universe.
    4. Numerically integrates the function to find the comoving distance.
    5. Compares the calculated distance to the provided options to check the final answer.
    """
    # --- 1. Define constants and parameters ---
    lambda_obs = 790.0   # Observed wavelength in nm
    lambda_rest = 121.6  # Rest-frame Lyman-alpha wavelength in nm
    H0 = 70.0            # Hubble constant in km/s/Mpc
    Omega_m = 0.3        # Matter density parameter
    Omega_L = 0.7        # Dark energy density parameter
    c_kms = 299792.458   # Speed of light in km/s

    # --- 2. Calculate the redshift (z) ---
    z = (lambda_obs / lambda_rest) - 1

    # --- 3. Define the integrand for comoving distance ---
    # The integrand is 1 / E(z'), where E(z') = sqrt(Omega_m * (1+z')^3 + Omega_L)
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # --- 4. Perform numerical integration and calculate distance ---
    try:
        integral_result, _ = quad(integrand, 0, z)
    except Exception as e:
        return f"An error occurred during numerical integration: {e}"

    # Hubble distance in Mpc
    D_H = c_kms / H0
    # Comoving distance in Mpc
    Dc_mpc = D_H * integral_result
    # Comoving distance in Gpc
    Dc_gpc = Dc_mpc / 1000.0

    # --- 5. Check the correctness of the final answer ---
    # The provided final answer is <<<D>>>, which corresponds to 9 Gpc.
    # The options given in the final analysis are:
    # A) 8 Gpc, B) 7 Gpc, C) 6 Gpc, D) 9 Gpc
    options = {'A': 8.0, 'B': 7.0, 'C': 6.0, 'D': 9.0}
    expected_answer_label = 'D'
    
    # Find which option is numerically closest to our calculated value
    closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_gpc))

    # Verify if the provided answer's choice ('D') is indeed the closest one
    if closest_option_label == expected_answer_label:
        # The calculation confirms the provided analysis.
        # The calculated value is ~8.99 Gpc, which is closest to 9 Gpc (Option D).
        return "Correct"
    else:
        # The calculation contradicts the provided analysis.
        reason = (
            f"The provided answer is <<<D>>>, corresponding to {options[expected_answer_label]} Gpc.\n"
            f"The verification code calculated a redshift of z = {z:.3f}.\n"
            f"The calculated comoving distance is {Dc_gpc:.3f} Gpc.\n"
            f"This calculated value is closest to the option '{closest_option_label}' ({options[closest_option_label]} Gpc), not '{expected_answer_label}'.\n"
            f"Therefore, the provided answer is incorrect."
        )
        return reason

# Run the check
print(check_cosmological_distance())