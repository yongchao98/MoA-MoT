import numpy as np
from scipy.integrate import quad

def check_cosmology_answer():
    """
    This function verifies the comoving distance calculation for the given quasar.
    It follows these steps:
    1. Defines the physical and cosmological parameters from the question.
    2. Calculates the redshift 'z' based on the Lyman-alpha line hypothesis.
    3. Numerically integrates the Friedmann equation to find the comoving distance.
    4. Compares the calculated distance to the provided options to check the answer's correctness.
    """
    # --- 1. Define constants and parameters ---
    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    # Physical constant: rest-frame wavelength of Lyman-alpha line
    lambda_rest = 121.567  # nm

    # Cosmological parameters from the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter
    c = 299792.458  # Speed of light in km/s

    # Options from the question
    options = {'A': 9.0, 'B': 7.0, 'C': 8.0, 'D': 6.0}
    # The final answer provided by the LLM
    llm_answer_label = 'A'

    # --- 2. Calculate the redshift (z) ---
    try:
        z = (lambda_obs / lambda_rest) - 1
    except Exception as e:
        return f"Error in redshift calculation: {e}"

    # --- 3. Calculate the comoving distance (Dc) ---
    # Define the function E(z') for the integrand's denominator
    def E(z_prime):
        return np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Define the integrand for the comoving distance integral
    def integrand(z_prime):
        return 1.0 / E(z_prime)

    # Perform the numerical integration from z=0 to the calculated redshift
    try:
        integral_value, _ = quad(integrand, 0, z)
    except Exception as e:
        return f"Error during numerical integration: {e}"

    # Calculate the Hubble distance in Mpc
    Dh_mpc = c / H0

    # Calculate the comoving distance in Mpc and then convert to Gpc
    Dc_mpc = Dh_mpc * integral_value
    Dc_gpc = Dc_mpc / 1000.0

    # --- 4. Verify the correctness of the answer ---
    # Find which option is mathematically closest to our calculated value
    closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_gpc))

    # Check if the LLM's answer matches the closest calculated option
    if llm_answer_label == closest_option_label:
        # The LLM's answer is correct.
        # Let's add a check for the reasoning. The reasoning states the value is ~9.0 Gpc.
        # A small tolerance (e.g., 2%) around the option value is reasonable.
        if abs(Dc_gpc - options[llm_answer_label]) / options[llm_answer_label] < 0.02:
            return "Correct"
        else:
            # This case is unlikely but would mean the closest option is correct by chance.
            return (f"Correct. The calculated distance is {Dc_gpc:.2f} Gpc, which is closest to option {llm_answer_label} "
                    f"({options[llm_answer_label]} Gpc), but the value is not as close as the reasoning implies.")
    else:
        # The LLM's answer is incorrect.
        return (f"Incorrect. The provided answer is {llm_answer_label} ({options[llm_answer_label]} Gpc), "
                f"but the calculated comoving distance is {Dc_gpc:.2f} Gpc. "
                f"This value is closest to option {closest_option_label} ({options[closest_option_label]} Gpc).")

# Run the check and print the result
print(check_cosmology_answer())