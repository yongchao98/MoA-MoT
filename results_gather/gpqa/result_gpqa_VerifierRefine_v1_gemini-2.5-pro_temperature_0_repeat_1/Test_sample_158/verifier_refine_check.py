import numpy as np
from scipy.integrate import quad
import sys
import io

def check_answer():
    """
    This function checks the correctness of the provided LLM answer for the cosmology problem.
    It recalculates the comoving distance based on the parameters given in the question
    and compares the result to the provided options.
    """
    # --- 1. Define constants and parameters from the question ---

    # Physical constants
    # Using a more precise value for the speed of light
    c_kms = 299792.458  # Speed of light in km/s
    # Observed wavelength from the problem statement
    lambda_obs_nm = 790.0  # nm
    # Rest-frame wavelength of the Lyman-alpha transition, a standard value in astronomy
    lambda_rest_lya_nm = 121.567 # nm

    # Cosmological parameters from the Lambda-CDM model in the question
    H0 = 70.0      # Hubble constant in km/s/Mpc
    Omega_m0 = 0.3   # Matter density parameter
    Omega_L0 = 0.7   # Dark energy density parameter

    # The provided answer to check
    llm_answer_option = 'C'
    options = {'A': 7.0, 'B': 6.0, 'C': 9.0, 'D': 8.0}


    # --- 2. Perform the step-by-step calculation ---

    # Step A: Verify the problem constraints
    # The problem states a flat universe. Let's check if the parameters are consistent.
    # For a flat universe, the total density parameter (Omega_total) should be 1.
    Omega_total = Omega_m0 + Omega_L0
    if not np.isclose(Omega_total, 1.0):
        return (f"Constraint check failed: The universe is stated to be flat, but the given density "
                f"parameters do not sum to 1 (Ω_m + Ω_Λ = {Omega_m0} + {Omega_L0} = {Omega_total}).")

    # Step B: Calculate the redshift (z)
    # The core assumption, as correctly identified in the LLM's analysis, is that the
    # peak at 790 nm followed by a flux drop is the redshifted Lyman-alpha line.
    # The formula for redshift is 1 + z = λ_obs / λ_rest
    try:
        z = (lambda_obs_nm / lambda_rest_lya_nm) - 1
    except ZeroDivisionError:
        return "Calculation Error: Rest-frame wavelength cannot be zero."

    # Step C: Define the function for numerical integration
    # The comoving distance D_c is calculated by integrating over redshift.
    # D_c = D_H * ∫[0 to z] dz' / E(z')
    # where D_H = c / H0 is the Hubble distance, and
    # E(z') = sqrt(Ω_m,₀ * (1+z')³ + Ω_Λ,₀) for a flat ΛCDM universe.

    def E(z_prime):
        """Dimensionless Hubble parameter for a flat Lambda-CDM universe."""
        return np.sqrt(Omega_m0 * (1 + z_prime)**3 + Omega_L0)

    def integrand(z_prime):
        """The function to be integrated: 1 / E(z')."""
        return 1.0 / E(z_prime)

    # Step D: Calculate the Hubble Distance (D_H)
    D_H_Mpc = c_kms / H0

    # Step E: Perform the numerical integration
    # We use scipy.integrate.quad for a precise numerical result.
    integral_result, integral_error = quad(integrand, 0, z)

    # Step F: Calculate the final comoving distance (D_c)
    Dc_Mpc = D_H_Mpc * integral_result
    Dc_Gpc = Dc_Mpc / 1000.0

    # --- 3. Compare the calculated result with the provided answer ---

    # Find which of the given options is closest to our calculated value.
    distances_to_options = {key: abs(val - Dc_Gpc) for key, val in options.items()}
    closest_option_key = min(distances_to_options, key=distances_to_options.get)

    # Check if the LLM's chosen option is the closest one.
    if closest_option_key == llm_answer_option:
        return "Correct"
    else:
        # If not, provide a detailed reason for the discrepancy.
        reason = (
            f"The provided answer '{llm_answer_option}' ({options[llm_answer_option]} Gpc) is incorrect.\n"
            f"The calculation based on the given parameters yields a comoving distance of {Dc_Gpc:.3f} Gpc.\n\n"
            f"Here is a summary of the calculation:\n"
            f"1. Redshift (z): Based on λ_obs={lambda_obs_nm} nm and λ_rest(Lyα)={lambda_rest_lya_nm} nm, the redshift is z = {z:.4f}.\n"
            f"2. Hubble Distance (D_H): D_H = c/H₀ = {c_kms:.1f} / {H0} ≈ {D_H_Mpc:.2f} Mpc.\n"
            f"3. Comoving Distance Integral: The integral of dz'/E(z') from 0 to {z:.4f} is ≈ {integral_result:.4f}.\n"
            f"4. Final Comoving Distance (D_c): D_c = D_H * integral ≈ {Dc_Mpc:.2f} Mpc = {Dc_Gpc:.3f} Gpc.\n\n"
            f"The calculated value of {Dc_Gpc:.3f} Gpc is closest to option '{closest_option_key}' ({options[closest_option_key]} Gpc), not '{llm_answer_option}'."
        )
        return reason

# Execute the check and print the result
# This captures the output of the function for display
# In a real application, you would just call `check_answer()`
# and use its return value.
output = check_answer()
print(output)
