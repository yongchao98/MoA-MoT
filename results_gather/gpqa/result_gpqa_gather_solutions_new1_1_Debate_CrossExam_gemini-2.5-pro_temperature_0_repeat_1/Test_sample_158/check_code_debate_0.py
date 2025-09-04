import numpy as np
from scipy.integrate import quad
import sys

def check_answer():
    """
    This function checks the correctness of the provided LLM answer for the cosmology question.
    It recalculates the comoving distance based on the problem statement and compares it
    to the LLM's reasoning and final choice.
    """
    # 1. Define the constants and parameters from the question
    # Physical constants
    c = 299792.458  # Speed of light in km/s

    # Cosmological parameters from the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km s^-1 Mpc^-1
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter

    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    
    # The problem states a "peak" with a flux drop at shorter wavelengths.
    # This is the classic signature of the Lyman-alpha (Ly-α) emission line.
    lambda_rest_lya = 121.567  # Rest-frame Ly-α wavelength in nm

    # The final answer provided by the LLM
    llm_answer = "D"
    options = {"A": 8.0, "B": 7.0, "C": 6.0, "D": 9.0}

    # 2. Step 1 from the LLM's reasoning: Calculate the redshift (z)
    try:
        z = (lambda_obs / lambda_rest_lya) - 1
    except Exception as e:
        return f"Error in redshift calculation: {e}"

    expected_z = 5.50
    if not np.isclose(z, expected_z, atol=0.01):
        return f"Redshift calculation is incorrect. Expected ~{expected_z}, but calculated {z:.4f}."

    # 3. Step 2 from the LLM's reasoning: Calculate the comoving distance (Dc)
    # The formula for comoving distance in a flat universe is:
    # Dc = (c / H0) * integral from 0 to z of [1 / E(z')] dz'
    # where E(z') = sqrt(Omega_m * (1 + z')^3 + Omega_L)

    def E(z_prime):
        return np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    def integrand(z_prime):
        return 1.0 / E(z_prime)

    try:
        # Perform the numerical integration
        integral_value, integral_error = quad(integrand, 0, z)
        
        # Calculate the Hubble Distance in Mpc
        Dh_mpc = c / H0
        
        # Calculate the comoving distance in Mpc and convert to Gpc
        Dc_mpc = Dh_mpc * integral_value
        Dc_gpc = Dc_mpc / 1000.0
    except Exception as e:
        return f"Error during comoving distance calculation: {e}"

    # The LLM's reasoning states the result is ~8.99 Gpc. Let's check.
    expected_Dc_gpc = 8.99
    if not np.isclose(Dc_gpc, expected_Dc_gpc, atol=0.05):
        return f"Comoving distance calculation is incorrect. The LLM claims ~{expected_Dc_gpc} Gpc, but the code calculated {Dc_gpc:.4f} Gpc."

    # 4. Step 3 from the LLM's reasoning: Compare with options
    # Find the option closest to our calculated value
    closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_gpc))

    if closest_option_label != llm_answer:
        return f"The final answer is incorrect. The calculated distance is {Dc_gpc:.4f} Gpc, which is closest to option {closest_option_label} ({options[closest_option_label]} Gpc), but the LLM chose option {llm_answer}."

    # 5. Final check: All steps match the LLM's reasoning and conclusion.
    return "Correct"

# Run the check
try:
    result = check_answer()
    print(result)
except ImportError:
    print("Could not run the check because required libraries (numpy, scipy) are not installed.")
except Exception as e:
    print(f"An unexpected error occurred: {e}")