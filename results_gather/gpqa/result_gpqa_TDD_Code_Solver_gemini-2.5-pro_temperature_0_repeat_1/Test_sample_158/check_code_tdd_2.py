import numpy as np
from scipy.integrate import quad

def check_answer():
    """
    Checks the correctness of the provided answer by recalculating the comoving distance
    based on the assumptions made in the answer.
    """
    # --- Define problem constraints and constants ---
    # Cosmological parameters from the question
    H0_km_s_Mpc = 70.0
    Omega_m0 = 0.3
    Omega_L0 = 0.7
    
    # Physical constants
    c_km_s = 299792.458  # Speed of light in km/s
    Mpc_per_Gpc = 1000.0

    # Observation from the question
    lambda_obs_nm = 790.0

    # --- Assumption from the provided answer ---
    # The answer assumes the spectral line is Carbon IV (C IV)
    lambda_rest_nm = 154.9
    
    # The chosen option in the answer
    chosen_option = 'A'
    options = {'A': 8.0, 'B': 9.0, 'C': 7.0, 'D': 6.0}

    # --- Step 1: Calculate Redshift (z) ---
    # z = (lambda_observed / lambda_rest) - 1
    if lambda_obs_nm < lambda_rest_nm:
        return "Calculation Error: Observed wavelength cannot be less than rest wavelength."
    z = (lambda_obs_nm / lambda_rest_nm) - 1.0

    # --- Step 2: Define the function to calculate comoving distance ---
    # For a flat Lambda-CDM universe, the comoving distance D_c is given by:
    # D_c = D_H * integral from 0 to z of dz' / E(z')
    # where D_H = c / H0 is the Hubble distance and E(z') = sqrt(Omega_m0*(1+z')^3 + Omega_L0)
    
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m0 * (1.0 + z_prime)**3 + Omega_L0)

    # --- Step 3: Perform the calculation ---
    try:
        hubble_distance_Mpc = c_km_s / H0_km_s_Mpc
        integral_part, integral_error = quad(integrand, 0, z)
        comoving_distance_Mpc = hubble_distance_Mpc * integral_part
        comoving_distance_Gpc = comoving_distance_Mpc / Mpc_per_Gpc
    except Exception as e:
        return f"An error occurred during integration: {e}"

    # --- Step 4: Verify the result against the answer ---
    # The answer's code calculates a distance of ~8.104 Gpc. Let's check our result.
    # A small tolerance is used for floating point comparison.
    if not np.isclose(comoving_distance_Gpc, 8.104, atol=0.01):
        return (f"Calculation mismatch: The provided answer's calculation leads to ~8.104 Gpc. "
                f"Our independent calculation under the same C IV assumption yields {comoving_distance_Gpc:.3f} Gpc. "
                f"This indicates a potential error in one of the calculations.")

    # Find the closest option to the calculated distance
    closest_option_key = min(options, key=lambda k: abs(options[k] - comoving_distance_Gpc))

    # Check if the closest option matches the answer's choice
    if closest_option_key == chosen_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculation based on the C IV assumption is correct and yields a comoving distance of "
                f"{comoving_distance_Gpc:.3f} Gpc. However, this value is closest to option {closest_option_key} ({options[closest_option_key]} Gpc), "
                f"not the chosen option {chosen_option} ({options[chosen_option]} Gpc).")

# Run the check
result = check_answer()
print(result)