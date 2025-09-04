import numpy as np
from scipy.integrate import quad

def check_answer():
    """
    This function checks the correctness of the provided answer by:
    1. Calculating the redshift from the spectral data.
    2. Calculating the comoving distance using the given cosmological model.
    3. Comparing the calculated distance to the options to see if the provided answer is the closest one.
    """
    # --- Define Constants and Parameters from the Question ---

    # Cosmological parameters
    H0 = 70.0  # Hubble constant in km s^-1 Mpc^-1
    Omega_m0 = 0.3  # Matter density parameter
    Omega_L0 = 0.7  # Dark energy density parameter (Lambda)
    # Constraint: Universe is flat, so Omega_k = 1 - Omega_m0 - Omega_L0 = 0. This is satisfied.

    # Physical constants and observational data
    c = 299792.458      # Speed of light in km/s
    lambda_obs = 790.0  # Observed wavelength in nm
    # Rest-frame wavelength of the Lyman-alpha transition
    lambda_rest_lya = 121.567 # in nm

    # Provided answer from the LLM
    llm_answer_option = 'D'
    options = {'A': 7.0, 'B': 6.0, 'C': 8.0, 'D': 9.0}
    llm_answer_value = options[llm_answer_option]

    # --- Step 1: Calculate the Redshift (z) ---
    # The spectral feature is identified as the redshifted Lyman-alpha line.
    # The formula for redshift is z = (λ_obs / λ_rest) - 1
    try:
        z = (lambda_obs / lambda_rest_lya) - 1
    except ZeroDivisionError:
        return "Error: Rest-frame wavelength cannot be zero."

    # --- Step 2: Calculate the Comoving Distance (Dc) ---
    # For a flat ΛCDM universe, the comoving distance is given by the integral:
    # Dc = (c/H0) * integral from 0 to z of [dz' / sqrt(Ω_m(1+z')^3 + Ω_Λ)]

    # Define the integrand E(z')^-1 where E(z') = H(z')/H0
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m0 * (1 + z_prime)**3 + Omega_L0)

    # Calculate the Hubble Distance in Mpc
    Dh_mpc = c / H0

    # Perform the numerical integration from z=0 to the calculated redshift
    try:
        integral_val, integral_err = quad(integrand, 0, z)
    except Exception as e:
        return f"Error during numerical integration: {e}"

    # Calculate the comoving distance in Mpc and then convert to Gpc
    Dc_mpc = Dh_mpc * integral_val
    Dc_gpc = Dc_mpc / 1000.0

    # --- Step 3: Check the Correctness of the LLM's Answer ---
    # Find which of the given options is closest to our calculated value.
    
    # Calculate the absolute difference between our result and each option
    differences = {key: abs(value - Dc_gpc) for key, value in options.items()}
    
    # Find the option with the minimum difference
    closest_option_key = min(differences, key=differences.get)

    # Check if the LLM's chosen option is the closest one
    if closest_option_key == llm_answer_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculated comoving distance is approximately {Dc_gpc:.2f} Gpc. "
                f"This value is closest to option {closest_option_key} ({options[closest_option_key]} Gpc), "
                f"not the provided answer {llm_answer_option} ({llm_answer_value} Gpc).")

# Run the check
result = check_answer()
print(result)