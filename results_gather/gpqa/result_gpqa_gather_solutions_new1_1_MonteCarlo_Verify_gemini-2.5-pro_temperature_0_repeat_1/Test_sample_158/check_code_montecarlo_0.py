import numpy as np
from scipy.integrate import quad

def check_cosmology_answer():
    """
    This function verifies the calculation for the comoving distance to a quasar.
    It follows the steps outlined in the provided analysis:
    1. Interpret the spectral feature as the Lyman-alpha line.
    2. Calculate the redshift (z).
    3. Numerically integrate to find the comoving distance (Dc) for the given Lambda-CDM model.
    4. Compare the calculated distance to the provided options to verify the final answer.
    """
    # --- Define constants and parameters from the question ---
    # Physical constants
    c = 299792.458  # Speed of light in km/s
    
    # Cosmological parameters
    H0 = 70.0       # Hubble constant in km/s/Mpc
    Omega_m = 0.3   # Matter density parameter
    Omega_L = 0.7   # Dark energy density parameter

    # Observational data and physical interpretation
    lambda_obs = 790.0      # Observed wavelength in nm
    # The analysis correctly identifies the feature as the Lyman-alpha (Lyα) line.
    lambda_rest = 121.567 # Rest-frame Lyα wavelength in nm

    # The final answer from the analysis is 'D', which corresponds to 9 Gpc.
    expected_answer_label = 'D'
    options = {'A': 6.0, 'B': 7.0, 'C': 8.0, 'D': 9.0}

    # --- Step 1: Calculate the redshift (z) ---
    try:
        z = (lambda_obs / lambda_rest) - 1
    except Exception as e:
        return f"Error during redshift calculation: {e}"

    # --- Step 2: Calculate the comoving distance (Dc) ---
    # Define the function E(z), the dimensionless Hubble parameter for a flat universe
    def E(z_prime):
        return np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Define the integrand for the comoving distance integral
    def integrand(z_prime):
        return 1.0 / E(z_prime)

    try:
        # Perform the numerical integration from 0 to z
        integral_value, _ = quad(integrand, 0, z)
    except Exception as e:
        return f"Error during numerical integration: {e}"

    # Calculate the Hubble distance in Mpc
    Dh_Mpc = c / H0

    # Calculate the comoving distance in Mpc and convert to Gpc
    Dc_Gpc = (Dh_Mpc * integral_value) / 1000.0

    # --- Step 3: Verify the final answer ---
    # Check if the primary interpretation (Lyα line) is consistent with the answer.
    if not (lambda_rest > 120 and lambda_rest < 122):
        return "Incorrect. The analysis is based on the Lyman-alpha line (λ_rest ≈ 121.6 nm), but the code is not using this value."

    # Check if the calculated redshift is consistent with the analysis.
    if not np.isclose(z, 5.5, atol=0.01):
        return f"Incorrect. The calculated redshift is z={z:.2f}, which deviates from the expected z≈5.5."

    # Check if the calculated distance is consistent with the analysis.
    if not np.isclose(Dc_Gpc, 8.98, atol=0.05):
        return f"Incorrect. The calculated comoving distance is {Dc_Gpc:.2f} Gpc, which deviates from the expected value of ~8.98 Gpc."

    # Find which option is mathematically closest to the calculated distance.
    closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_Gpc))

    # Check if the closest option matches the provided answer.
    if closest_option_label == expected_answer_label:
        return "Correct"
    else:
        return (f"Incorrect. The calculated comoving distance is {Dc_Gpc:.2f} Gpc. "
                f"The closest option is {closest_option_label} ({options[closest_option_label]} Gpc), "
                f"but the provided answer was {expected_answer_label} ({options[expected_answer_label]} Gpc).")

# Run the check
result = check_cosmology_answer()
print(result)