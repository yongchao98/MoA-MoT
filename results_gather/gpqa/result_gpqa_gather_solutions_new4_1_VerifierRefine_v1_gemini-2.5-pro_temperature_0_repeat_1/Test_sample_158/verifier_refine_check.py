import numpy as np
from scipy.integrate import quad

def check_correctness():
    """
    This function checks the correctness of the provided answer to the cosmology question.
    It recalculates the comoving distance based on the problem's parameters and compares
    the result to the provided answer.
    """
    try:
        # --- Step 1: Define Constants and Parameters ---
        # Observational data from the question
        lambda_obs = 790.0  # Observed wavelength in nm
        # Use a precise value for the rest-frame Lyman-alpha wavelength
        lambda_rest_lya = 121.567  # in nm

        # Cosmological parameters from the question
        H0 = 70.0  # Hubble constant in km/s/Mpc
        Omega_m = 0.3  # Matter density parameter
        Omega_L = 0.7  # Dark energy density parameter
        
        # Physical constants
        c = 299792.458  # Speed of light in km/s

        # --- Step 2: Calculate the Redshift (z) ---
        # The redshift is calculated from the observed and rest wavelengths
        z = (lambda_obs / lambda_rest_lya) - 1

        # --- Step 3: Calculate the Comoving Distance (D_c) ---
        # The comoving distance requires a numerical integration.
        # First, define the dimensionless Hubble parameter E(z) for a flat Lambda-CDM universe.
        def E(z_prime):
            return np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

        # The integrand for the comoving distance integral is 1/E(z).
        def integrand(z_prime):
            return 1.0 / E(z_prime)

        # Perform the numerical integration from z=0 to the calculated redshift z.
        # quad returns a tuple (integral_value, estimated_error)
        integral_value, _ = quad(integrand, 0, z)

        # Calculate the Hubble distance in Mpc
        D_H = c / H0

        # Calculate the comoving distance in Mpc
        D_C_Mpc = D_H * integral_value
        
        # Convert the final result to Gigaparsecs (Gpc)
        D_C_Gpc = D_C_Mpc / 1000.0

        # --- Step 4: Compare the Calculated Result with the Provided Answer ---
        # The multiple-choice options from the question
        options = {'A': 8.0, 'B': 7.0, 'C': 9.0, 'D': 6.0}
        
        # The final answer provided in the prompt to be checked is 'C'
        provided_answer_label = 'C'
        
        # Find which option is mathematically closest to our calculated value
        distances = {label: abs(D_C_Gpc - value) for label, value in options.items()}
        closest_option_label = min(distances, key=distances.get)

        # Check if the closest option matches the provided answer
        if closest_option_label == provided_answer_label:
            return "Correct"
        else:
            reason = (
                f"Incorrect.\n"
                f"The provided answer is '{provided_answer_label}' ({options[provided_answer_label]} Gpc), but the calculation shows this is not the closest option.\n"
                f"Calculation Steps:\n"
                f"1. The redshift (z) is calculated as ({lambda_obs} nm / {lambda_rest_lya} nm) - 1 = {z:.4f}.\n"
                f"2. The comoving distance is calculated by numerically integrating the Friedmann equation, yielding {D_C_Gpc:.4f} Gpc.\n"
                f"3. This calculated value of {D_C_Gpc:.4f} Gpc is closest to option '{closest_option_label}' ({options[closest_option_label]} Gpc).\n"
                f"Therefore, the provided answer '{provided_answer_label}' is incorrect."
            )
            return reason

    except ImportError:
        return "Could not run the check because required libraries (numpy, scipy) are not installed."
    except Exception as e:
        return f"An error occurred during the check: {e}"

# Execute the check and print the result
result = check_correctness()
print(result)