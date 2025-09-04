import numpy as np
from scipy.integrate import quad

def check_quasar_distance():
    """
    This function checks the correctness of the final answer by performing the
    necessary calculations from first principles.
    """
    # --- Step 1: Define constants and parameters from the question ---
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest = 121.567  # Rest-frame Lyman-alpha wavelength in nm
    
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter
    
    c = 299792.458  # Speed of light in km/s

    # --- Step 2: Calculate the redshift (z) ---
    # z = (lambda_obs / lambda_rest) - 1
    z = (lambda_obs / lambda_rest) - 1
    
    # --- Step 3: Calculate the comoving distance ---
    # The integrand for the comoving distance calculation in a flat Lambda-CDM universe
    def integrand(z_prime):
        E_z = np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)
        return 1.0 / E_z

    # Perform the numerical integration from z=0 to the calculated redshift z
    integral_val, _ = quad(integrand, 0, z)
    
    # Calculate the Hubble distance in Mpc (D_H = c / H0)
    hubble_distance_mpc = c / H0
    
    # Calculate the comoving distance in Mpc
    comoving_distance_mpc = hubble_distance_mpc * integral_val
    
    # Convert the result to Gigaparsecs (Gpc)
    calculated_distance_gpc = comoving_distance_mpc / 1000.0

    # --- Step 4: Check the correctness of the provided final answer ---
    # The final answer provided is <<<A>>>, which corresponds to 9 Gpc in its option list.
    # Option list from the final answer's reasoning: A) 9 Gpc, B) 7 Gpc, C) 6 Gpc, D) 8 Gpc
    options = {'A': 9.0, 'B': 7.0, 'C': 6.0, 'D': 8.0}
    final_answer_choice = 'A'
    final_answer_value = options[final_answer_choice]

    # Find which of the given options is closest to our calculated value
    option_values = np.array(list(options.values()))
    closest_option_value = option_values[np.argmin(np.abs(option_values - calculated_distance_gpc))]

    # The calculation is correct if the provided answer's value is the closest one.
    # We use a tolerance for floating point comparisons.
    if not np.isclose(final_answer_value, closest_option_value):
        return (f"Incorrect. The provided answer is {final_answer_value} Gpc. "
                f"However, the calculated comoving distance is approximately {calculated_distance_gpc:.2f} Gpc. "
                f"The closest option to the calculated value is {closest_option_value} Gpc.")

    # A final check to ensure the calculation itself is sound. The discrepancy in the
    # candidate answers was between ~8 Gpc and ~9 Gpc. Our calculation should be close to 9 Gpc.
    if np.isclose(calculated_distance_gpc, 9.0, atol=0.2):
        return "Correct"
    else:
        # This case would trigger if the provided answer was 9 Gpc, but our calculation
        # yielded a significantly different result (e.g., 8 Gpc).
        return (f"Incorrect. The provided answer is {final_answer_value} Gpc, but the verification "
                f"calculation yields a comoving distance of {calculated_distance_gpc:.2f} Gpc. "
                f"This contradicts the provided answer.")

# Run the check
result = check_quasar_distance()
print(result)