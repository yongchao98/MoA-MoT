import numpy as np
from scipy.integrate import quad

def check_correctness_of_cosmology_answer():
    """
    This function verifies the calculation of the comoving distance to a quasar
    based on its observed spectral features and a given cosmological model.

    It follows these steps:
    1.  Calculates the redshift (z) from the observed Lyman-alpha line.
    2.  Numerically integrates the Friedmann equation for a flat Lambda-CDM universe
        to find the comoving distance.
    3.  Compares the calculated distance to the given options to verify the provided answer.
    """
    # --- Problem Parameters ---
    # Cosmological model parameters
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter

    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest_lya = 121.567  # Rest-frame Lyman-alpha wavelength in nm

    # Physical constants
    c = 299792.458  # Speed of light in km/s

    # --- Step 1: Calculate the redshift (z) ---
    # The redshift is determined by the stretching of the Lyman-alpha line's wavelength.
    try:
        z = (lambda_obs / lambda_rest_lya) - 1
    except Exception as e:
        return f"Error during redshift calculation: {e}"

    # --- Step 2: Calculate the comoving distance ---
    # The comoving distance in a flat Lambda-CDM universe is given by the integral:
    # D_c = (c / H0) * integral from 0 to z of [dz' / E(z')]
    # where E(z) = sqrt(Omega_m * (1+z)^3 + Omega_L)

    # Define the function E(z), the dimensionless Hubble parameter
    def E(z_prime):
        return np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Define the integrand for the comoving distance integral
    def integrand(z_prime):
        return 1.0 / E(z_prime)

    # Perform the numerical integration from 0 to z
    try:
        integral_value, _ = quad(integrand, 0, z)
    except Exception as e:
        return f"Error during numerical integration: {e}"

    # Calculate the Hubble distance in Mpc
    hubble_distance_mpc = c / H0

    # Calculate the comoving distance in Mpc and then convert to Gpc
    comoving_distance_mpc = hubble_distance_mpc * integral_value
    comoving_distance_gpc = comoving_distance_mpc / 1000.0

    # --- Step 3: Verify the final answer ---
    # The provided final answer is 'D', which corresponds to 9 Gpc from the options.
    # Options from the question: A) 6 Gpc, B) 8 Gpc, C) 7 Gpc, D) 9 Gpc
    options = {'A': 6, 'B': 8, 'C': 7, 'D': 9}
    final_answer_label = 'D'
    expected_value = options[final_answer_label]

    # Find which option value is mathematically closest to our calculated distance
    closest_option_value = min(options.values(), key=lambda val: abs(val - comoving_distance_gpc))

    # Check if the expected answer value (from 'D') is the closest one
    if np.isclose(closest_option_value, expected_value, atol=0.1):
        return "Correct"
    else:
        return (f"Incorrect. The calculated comoving distance is approximately {comoving_distance_gpc:.2f} Gpc. "
                f"This value is closest to the option {closest_option_value} Gpc. "
                f"The provided answer was '{final_answer_label}', which corresponds to {expected_value} Gpc, but this is not the closest option.")

# Run the check and print the result
result = check_correctness_of_cosmology_answer()
print(result)