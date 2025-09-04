import numpy as np
from scipy.integrate import quad

def check_quasar_distance():
    """
    This function checks the correctness of the LLM's answer by recalculating the comoving distance.
    It follows these steps:
    1.  Define the cosmological and physical constants from the problem statement.
    2.  Calculate the redshift (z) based on the observed Lyman-alpha line.
    3.  Calculate the theoretical comoving distance using the provided Lambda-CDM model parameters.
    4.  Compare the calculated distance with the provided answer (Option A: 9 Gpc).
    """
    # --- Step 1: Define constants and parameters ---
    # Cosmological parameters from the question
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_Lambda = 0.7  # Dark energy density parameter

    # Physical constants
    c = 299792.458  # Speed of light in km/s

    # Observational data from the question
    lambda_obs = 790.0  # Observed wavelength in nm
    # The spectral feature is identified as the Lyman-alpha (Lyα) line
    lambda_rest_lya = 121.567  # Rest-frame wavelength of Lyα in nm

    # The answer to check
    llm_answer_option = 'A'
    options = {'A': 9.0, 'B': 7.0, 'C': 8.0, 'D': 6.0} # in Gpc
    llm_answer_value = options[llm_answer_option]

    # --- Step 2: Calculate the redshift (z) ---
    # The redshift is calculated from the shift in the spectral line.
    # 1 + z = λ_obs / λ_rest
    z = (lambda_obs / lambda_rest_lya) - 1

    # --- Step 3: Calculate the comoving distance ---
    # The formula for comoving distance (D_c) in a flat Lambda-CDM universe is:
    # D_c = D_H * ∫[0 to z] dz' / E(z')
    # where D_H = c / H₀ is the Hubble distance, and
    # E(z') = sqrt(Ω_m * (1 + z')³ + Ω_Λ)

    # Define the integrand, 1 / E(z')
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_Lambda)

    # Calculate the Hubble Distance in Mpc
    hubble_distance_mpc = c / H0

    # Perform the numerical integration from z=0 to the calculated redshift z
    try:
        integral_value, _ = quad(integrand, 0, z)
    except Exception as e:
        return f"An error occurred during the numerical integration step: {e}"

    # Calculate the comoving distance in Mpc and then convert to Gpc
    comoving_distance_mpc = hubble_distance_mpc * integral_value
    calculated_distance_gpc = comoving_distance_mpc / 1000.0

    # --- Step 4: Compare calculated distance with the LLM's answer ---
    # Find the option that is numerically closest to our calculation
    closest_option = min(options, key=lambda k: abs(options[k] - calculated_distance_gpc))

    # Check if the LLM's answer is indeed the closest option
    if llm_answer_option != closest_option:
        return (f"Incorrect. The calculated comoving distance is {calculated_distance_gpc:.3f} Gpc. "
                f"This is closest to option {closest_option} ({options[closest_option]} Gpc), "
                f"not the provided answer {llm_answer_option} ({llm_answer_value} Gpc).")

    # Check if the value is reasonably close (e.g., within 0.2 Gpc)
    difference = abs(llm_answer_value - calculated_distance_gpc)
    if difference > 0.2:
        return (f"Incorrect. Although option {llm_answer_option} is the closest, the calculated distance "
                f"({calculated_distance_gpc:.3f} Gpc) differs from the answer's value ({llm_answer_value} Gpc) "
                f"by {difference:.3f} Gpc, which is a significant margin.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check and print the result.
result = check_quasar_distance()
print(result)