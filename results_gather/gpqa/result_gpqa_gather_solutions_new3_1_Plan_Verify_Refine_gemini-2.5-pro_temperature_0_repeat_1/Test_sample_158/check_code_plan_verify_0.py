import numpy as np
from scipy.integrate import quad

def check_correctness():
    """
    This function checks the correctness of the provided answer by recalculating the physical quantities.
    
    The steps are:
    1.  Identify the spectral feature as the Lyman-alpha (Lyα) line and calculate the redshift (z).
    2.  Define the cosmological parameters from the question.
    3.  Numerically integrate the Friedmann equation to find the comoving distance (Dc).
    4.  Compare the calculated distance to the given options and check if it matches the provided answer.
    """
    
    # --- Step 1: Define constants and calculate redshift ---
    
    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    
    # Physical constants and interpretation
    # The spectral feature is interpreted as the Lyman-alpha (Lyα) line.
    lambda_rest_lya = 121.567  # Rest-frame wavelength of Lyα in nm
    
    # Calculate redshift
    try:
        z = (lambda_obs / lambda_rest_lya) - 1
    except Exception as e:
        return f"Failed during redshift calculation: {e}"

    # --- Step 2: Define cosmological model and calculate distance ---
    
    # Cosmological parameters from the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km s^-1 Mpc^-1
    Omega_m = 0.3  # Matter density parameter
    Omega_lambda = 0.7  # Dark energy density parameter
    c = 299792.458  # Speed of light in km/s

    # Define the integrand for the comoving distance calculation.
    # This is the inverse of the dimensionless Hubble parameter E(z).
    def integrand(z_prime):
        E_z = np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)
        return 1.0 / E_z

    try:
        # Calculate the Hubble Distance in Mpc.
        hubble_distance_mpc = c / H0
        
        # Perform the numerical integration from z=0 to the calculated redshift.
        integral_val, _ = quad(integrand, 0, z)
        
        # Calculate the comoving distance in Mpc and convert to Gpc.
        dc_mpc = hubble_distance_mpc * integral_val
        dc_gpc = dc_mpc / 1000.0
    except Exception as e:
        return f"Failed during comoving distance calculation: {e}"

    # --- Step 3: Compare with options and the provided answer ---
    
    # Options from the question are: A) 9 Gpc, B) 7 Gpc, C) 8 Gpc, D) 6 Gpc
    options = {
        "A": 9.0,
        "B": 7.0,
        "C": 8.0,
        "D": 6.0
    }
    
    # The final answer provided by the LLM is 'A'.
    provided_answer_label = 'A'
    
    # Find the option that is numerically closest to our calculated distance.
    closest_option_label = min(options, key=lambda k: abs(options[k] - dc_gpc))

    # --- Step 4: Final check and return result ---
    
    # The primary assumption is that the 790 nm peak is the redshifted Lyα line.
    # This leads to a redshift z ≈ 5.5.
    # The comoving distance for this redshift in the given cosmology is calculated.
    # The final answer states the result is closest to 9 Gpc (Option A).
    # We check if our calculation confirms this.
    
    if closest_option_label == provided_answer_label:
        return "Correct"
    else:
        calculated_value = options[closest_option_label]
        reason = (
            f"The provided answer is '{provided_answer_label}' ({options[provided_answer_label]} Gpc). "
            f"However, the calculation based on the problem statement yields a comoving distance of approximately {dc_gpc:.2f} Gpc. "
            f"This calculated value is closest to option '{closest_option_label}' ({calculated_value} Gpc). "
            f"Therefore, the provided answer is incorrect."
        )
        return reason

# Execute the check and print the result
result = check_correctness()
print(result)