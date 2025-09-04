import numpy as np
from scipy.integrate import quad

def check_quasar_distance_answer():
    """
    Checks the correctness of the LLM's answer by calculating the comoving distance
    under different physical assumptions for the observed spectral line.
    """
    # --- Problem Constants and Cosmological Parameters ---
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m0 = 0.3  # Matter density parameter
    Omega_L0 = 0.7  # Dark energy density parameter
    # The universe is flat, so Omega_k = 1 - Omega_m0 - Omega_L0 = 0.

    c_km_s = 299792.458  # Speed of light in km/s
    Mpc_per_Gpc = 1000.0

    lambda_obs_nm = 790.0  # Observed wavelength in nm

    # --- Rest-frame wavelengths for possible spectral lines ---
    # LLM's assumption: Carbon IV
    lambda_rest_civ_nm = 154.9
    # Standard assumption: Lyman-alpha
    lambda_rest_lya_nm = 121.567

    # --- Function to calculate comoving distance ---
    def calculate_comoving_distance(z):
        """Calculates comoving distance in Gpc for a given redshift z."""
        if z <= 0:
            return 0.0
        
        # Integrand for comoving distance in a flat Lambda-CDM universe
        def integrand(z_prime):
            return 1.0 / np.sqrt(Omega_m0 * (1.0 + z_prime)**3 + Omega_L0)

        hubble_distance_Mpc = c_km_s / H0
        integral_part, _ = quad(integrand, 0, z)
        comoving_distance_Mpc = hubble_distance_Mpc * integral_part
        return comoving_distance_Mpc / Mpc_per_Gpc

    # --- Scenario 1: LLM's Assumption (C IV line) ---
    z_civ = (lambda_obs_nm / lambda_rest_civ_nm) - 1
    dist_civ_gpc = calculate_comoving_distance(z_civ)

    # --- Scenario 2: Standard Assumption (Lyman-alpha line) ---
    z_lya = (lambda_obs_nm / lambda_rest_lya_nm) - 1
    dist_lya_gpc = calculate_comoving_distance(z_lya)

    # --- Analysis ---
    llm_answer_option = 'C'
    options = {'A': 8, 'B': 9, 'C': 7, 'D': 6}
    llm_answer_value = options[llm_answer_option]

    # Check which option is closest for each scenario
    closest_option_civ = min(options, key=lambda k: abs(options[k] - dist_civ_gpc))
    closest_option_lya = min(options, key=lambda k: abs(options[k] - dist_lya_gpc))

    # The LLM's calculation is internally consistent.
    # It assumes C IV, calculates z ~ 4.1, gets D_c ~ 7.24 Gpc, and chooses C (7 Gpc).
    # Let's verify this logic chain.
    if closest_option_civ != llm_answer_option:
        return (f"Incorrect. The LLM's own logic is flawed. "
                f"Assuming the C IV line (z={z_civ:.2f}), the calculated distance is {dist_civ_gpc:.2f} Gpc. "
                f"The closest option should be {closest_option_civ}, but the LLM chose {llm_answer_option}.")

    # Now, evaluate the physical premise.
    # The description of the spectral feature strongly implies a Lyman-alpha break.
    # This is the standard interpretation in astrophysics.
    if closest_option_lya != llm_answer_option:
        return (f"Incorrect. The answer is based on a questionable physical assumption. "
                f"The problem describes a 'peak... and at shorter wavelengths the flux drops significantly', which is the classic signature of the Lyman-alpha break in a high-redshift quasar. "
                f"This implies the observed line is Lyman-alpha (rest wavelength ~121.6 nm), not Carbon IV.\n\n"
                f"1. LLM's Assumption (C IV line at 154.9 nm):\n"
                f"   - Redshift z = (790 / 154.9) - 1 = {z_civ:.3f}\n"
                f"   - Comoving Distance = {dist_civ_gpc:.3f} Gpc.\n"
                f"   - This is closest to option C (7 Gpc), which matches the LLM's answer.\n\n"
                f"2. Standard Physical Assumption (Lyman-alpha line at 121.6 nm):\n"
                f"   - Redshift z = (790 / 121.6) - 1 = {z_lya:.3f}\n"
                f"   - Comoving Distance = {dist_lya_gpc:.3f} Gpc.\n"
                f"   - This is closest to option A (8 Gpc).\n\n"
                f"Conclusion: While the LLM's calculation is correct for its chosen assumption, the assumption itself is physically less plausible than the standard interpretation, which leads to a different answer (A). Therefore, the final answer 'C' is incorrect.")

    return "Correct"

# Run the check
result = check_quasar_distance_answer()
print(result)