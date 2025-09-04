import numpy as np
from scipy.integrate import quad

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the comoving distance.
    """
    # --- Define Constants and Parameters from the Question ---

    # Physical constants
    c = 299792.458  # Speed of light in km/s
    
    # Cosmological parameters from the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_lambda = 0.7  # Dark energy density parameter

    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    
    # The spectral feature "peak ... with flux drop at shorter wavelengths" is the Lyman-alpha line.
    lambda_rest = 121.567 # Rest-frame wavelength of Lyman-alpha in nm

    # Multiple choice options from the question
    options = {
        'A': 6.0,
        'B': 7.0,
        'C': 9.0,
        'D': 8.0
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_key = 'C'

    # --- Step 1: Calculate the Redshift (z) ---
    try:
        z = (lambda_obs / lambda_rest) - 1
    except Exception as e:
        return f"Error during redshift calculation: {e}"

    # --- Step 2: Calculate the Comoving Distance ---
    
    # Define the integrand for the comoving distance calculation.
    # This is the inverse of the dimensionless Hubble parameter E(z).
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)

    try:
        # Calculate the Hubble Distance in Mpc
        hubble_distance_mpc = c / H0
        
        # Perform the numerical integration from z=0 to the calculated redshift
        integral_value, integral_error = quad(integrand, 0, z)
        
        # Calculate the comoving distance in Mpc and convert to Gpc
        comoving_distance_mpc = hubble_distance_mpc * integral_value
        comoving_distance_gpc = comoving_distance_mpc / 1000.0
    except Exception as e:
        return f"Error during comoving distance calculation: {e}"

    # --- Step 3: Compare the result with the options ---
    
    # Find the option closest to the calculated value
    closest_option_key = min(options, key=lambda k: abs(options[k] - comoving_distance_gpc))
    
    # --- Step 4: Validate the LLM's answer ---
    
    if closest_option_key == llm_answer_key:
        return "Correct"
    else:
        return (f"Incorrect. The calculation based on the problem's constraints yields a comoving distance "
                f"of {comoving_distance_gpc:.2f} Gpc. This value is closest to option "
                f"{closest_option_key} ({options[closest_option_key]} Gpc), but the provided answer was "
                f"{llm_answer_key} ({options[llm_answer_key]} Gpc).")

# Run the check and print the result
result = check_answer()
print(result)