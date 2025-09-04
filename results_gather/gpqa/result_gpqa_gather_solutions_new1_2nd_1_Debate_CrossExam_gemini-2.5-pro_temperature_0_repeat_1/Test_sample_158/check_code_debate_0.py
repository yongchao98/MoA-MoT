import numpy as np
from scipy import integrate
from scipy.constants import c

def check_cosmological_distance():
    """
    This function checks the correctness of the final answer by recalculating
    the comoving distance from the given observational data and cosmological model.
    """
    # --- Step 1: Define constants and parameters from the question ---
    lambda_obs_nm = 790.0      # Observed wavelength in nm
    lambda_rest_nm = 121.567   # Precise rest-frame Lyman-alpha wavelength in nm
    H0_km_s_Mpc = 70.0         # Hubble constant in km/s/Mpc
    Omega_m = 0.3              # Matter density parameter
    Omega_L = 0.7              # Dark energy density parameter (for a flat universe)

    # Physical constants
    c_km_s = c / 1000.0        # Speed of light in km/s

    # --- Step 2: Calculate the redshift (z) ---
    # z = (lambda_observed / lambda_rest) - 1
    redshift = (lambda_obs_nm / lambda_rest_nm) - 1

    # --- Step 3: Calculate the comoving distance ---
    # The comoving distance Dc is given by the integral of c/(H0 * E(z))
    # where E(z) = sqrt(Omega_m * (1+z)^3 + Omega_L) for a flat universe.
    
    # Define the function to be integrated: 1 / E(z')
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Perform the numerical integration from z'=0 to z'=redshift
    integral_val, _ = integrate.quad(integrand, 0, redshift)

    # Calculate the Hubble distance in Mpc: D_H = c / H0
    D_H_Mpc = c_km_s / H0_km_s_Mpc

    # Calculate the comoving distance in Mpc: Dc = D_H * integral
    Dc_Mpc = D_H_Mpc * integral_val

    # Convert the result to Gigaparsecs (Gpc)
    Dc_Gpc = Dc_Mpc / 1000.0

    # --- Step 4: Check the correctness of the provided answer ---
    # The provided answer is <<<A>>>, which corresponds to 9 Gpc.
    expected_answer_Gpc = 9.0
    
    # Check if the calculated value is close to the expected answer.
    # A tolerance of 2% (0.18 Gpc) is reasonable for this type of problem.
    if abs(Dc_Gpc - expected_answer_Gpc) / expected_answer_Gpc < 0.02:
        return "Correct"
    else:
        # If incorrect, provide a detailed reason.
        reason = (
            f"The provided answer 'A' corresponds to 9 Gpc, but the calculation yields a different result.\n"
            f"1. Based on the Lyα line (λ_rest={lambda_rest_nm} nm) and the observed peak (λ_obs={lambda_obs_nm} nm), the calculated redshift is z = {redshift:.3f}.\n"
            f"2. Using the given cosmological parameters (H₀={H0_km_s_Mpc}, Ω_m={Omega_m}, Ω_Λ={Omega_L}), the comoving distance is calculated.\n"
            f"3. The calculated comoving distance is {Dc_Gpc:.2f} Gpc.\n"
            f"4. This value ({Dc_Gpc:.2f} Gpc) is not consistent with the provided answer of {expected_answer_Gpc} Gpc."
        )
        # Determine which option is the closest to the calculated value.
        options = {'A': 9, 'B': 8, 'C': 6, 'D': 7}
        closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_Gpc))
        closest_option_value = options[closest_option_label]
        reason += f"\nThe calculated value is closest to option {closest_option_label} ({closest_option_value} Gpc)."
        return reason

# Execute the check and print the result.
result = check_cosmological_distance()
print(result)