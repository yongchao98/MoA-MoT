import numpy as np
from scipy.integrate import quad
from scipy.constants import c

def check_quasar_distance():
    """
    This function checks the correctness of the LLM's answer regarding the comoving distance to a quasar.
    It recalculates the redshift and the comoving distance based on the provided cosmological parameters.
    """
    # --- Parameters from the question ---
    lambda_obs_nm = 790.0  # Observed wavelength in nm
    lambda_rest_lya_nm = 121.6  # Rest-frame Lyman-alpha wavelength in nm
    H0_kms_Mpc = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_lambda = 0.7  # Dark energy density parameter
    
    # The LLM's chosen answer and calculated value
    llm_answer_option = 'B'
    llm_calculated_distance_gpc = 8.02
    
    # --- Step 1: Recalculate the Redshift (z) ---
    # The LLM correctly identifies the spectral feature as the redshifted Lyman-alpha line.
    # The redshift z is given by z = (lambda_obs / lambda_rest) - 1
    try:
        z = (lambda_obs_nm / lambda_rest_lya_nm) - 1
    except ZeroDivisionError:
        return "Error: Rest-frame wavelength cannot be zero."

    # Check if the LLM's redshift calculation is correct
    llm_z = 5.5
    if not np.isclose(z, llm_z, atol=0.01):
        return f"Incorrect redshift calculation. The LLM calculated z ~ {llm_z}, but the correct value is z = ({lambda_obs_nm} / {lambda_rest_lya_nm}) - 1 = {z:.4f}."

    # --- Step 2: Recalculate the Comoving Distance (Dc) ---
    # For a flat Lambda-CDM universe, the comoving distance is given by:
    # Dc = DH * integral from 0 to z of dz' / E(z')
    # where DH = c / H0 is the Hubble distance and E(z') = sqrt(Omega_m * (1+z')^3 + Omega_lambda)
    
    # Convert speed of light to km/s to match H0 units
    c_kms = c / 1000.0
    
    # Hubble distance in Mpc
    hubble_distance_mpc = c_kms / H0_kms_Mpc
    
    # Define the integrand for the comoving distance calculation
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)
        
    # Perform the numerical integration from 0 to z
    try:
        integral_result, integral_error = quad(integrand, 0, z)
    except Exception as e:
        return f"An error occurred during numerical integration: {e}"

    # Calculate the comoving distance in Mpc
    comoving_distance_mpc = hubble_distance_mpc * integral_result
    
    # Convert to Gpc
    comoving_distance_gpc = comoving_distance_mpc / 1000.0
    
    # --- Step 3: Compare with the LLM's answer ---
    # The LLM's answer is B, which corresponds to 8 Gpc.
    # The LLM's detailed calculation gives 8.02 Gpc.
    
    # Check if our calculated distance is close to the LLM's calculation and the chosen option.
    # We use a tolerance, e.g., 2%, to account for rounding in parameters or numerical precision.
    target_distance_gpc = 8.0
    tolerance = 0.1 # Gpc, a reasonable margin for this kind of problem
    
    if abs(comoving_distance_gpc - target_distance_gpc) <= tolerance:
        # The calculation is consistent with the LLM's result and the chosen option.
        return "Correct"
    else:
        # The calculation does not match.
        return (f"Incorrect. The final comoving distance calculation is wrong.\n"
                f"Given z = {z:.2f}, H0 = {H0_kms_Mpc}, Omega_m = {Omega_m}, and Omega_lambda = {Omega_lambda}, "
                f"the calculated comoving distance is {comoving_distance_gpc:.2f} Gpc.\n"
                f"The LLM's answer was {llm_answer_option} ({target_distance_gpc} Gpc), which is not consistent with the calculated value.")

# Run the check
result = check_quasar_distance()
print(result)