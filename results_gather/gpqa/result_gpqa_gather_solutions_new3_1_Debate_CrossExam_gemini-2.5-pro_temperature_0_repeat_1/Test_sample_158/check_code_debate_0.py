import numpy as np
from scipy.integrate import quad

def check_quasar_distance():
    """
    This function checks the correctness of the provided answer for the quasar distance problem.
    It performs the following steps:
    1.  Calculates the redshift (z) based on the Lyman-alpha (Lyα) line hypothesis.
    2.  Defines the cosmological parameters from the question.
    3.  Numerically integrates the Friedmann equation to find the comoving distance.
    4.  Compares the calculated distance to the given options to find the closest match.
    5.  Checks if this match corresponds to the provided final answer.
    """
    # --- Step 1: Define constants and parameters from the question ---
    # Physical constants and observational data
    c = 299792.458  # Speed of light in km/s
    lambda_obs = 790.0  # Observed wavelength in nm
    # The spectral feature is the Lyman-alpha (Lyα) line
    lambda_rest_lya = 121.567 # Rest-frame wavelength of Lyα in nm

    # Cosmological parameters for the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3
    Omega_lambda = 0.7

    # The final answer from the LLM to be checked
    llm_final_answer_option = 'D'
    options = {'A': 6.0, 'B': 9.0, 'C': 7.0, 'D': 8.0}

    # --- Step 2: Calculate redshift (z) ---
    # The most plausible interpretation is that the 790 nm peak is the redshifted Lyα line.
    z = (lambda_obs / lambda_rest_lya) - 1

    # --- Step 3: Calculate the comoving distance ---
    # Define the integrand for the comoving distance calculation.
    # This is 1/E(z), where E(z) is the dimensionless Hubble parameter.
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_lambda)

    # Calculate the Hubble distance in Mpc
    hubble_distance_mpc = c / H0

    # Perform the numerical integration from z=0 to the calculated redshift
    try:
        integral_value, integral_error = quad(integrand, 0, z)
    except Exception as e:
        return f"An error occurred during numerical integration: {e}"

    # Calculate the comoving distance in Mpc and convert to Gpc
    comoving_distance_mpc = hubble_distance_mpc * integral_value
    comoving_distance_gpc = comoving_distance_mpc / 1000.0

    # --- Step 4: Verify the final answer ---
    # Find which of the multiple-choice options is closest to our calculated value
    option_values = np.array(list(options.values()))
    option_labels = list(options.keys())
    closest_option_index = np.argmin(np.abs(option_values - comoving_distance_gpc))
    closest_option_label = option_labels[closest_option_index]

    # Check if the LLM's final answer matches the result of our calculation
    if llm_final_answer_option == closest_option_label:
        return "Correct"
    else:
        reason = (
            f"The final answer is incorrect. The analysis provided in the prompt is correct, but the final letter chosen is wrong in some of the candidate answers.\n"
            f"1. The spectral feature at 790 nm is the Lyman-alpha line (rest wavelength ~121.6 nm), which gives a redshift z ≈ {z:.2f}.\n"
            f"2. Using the given cosmological parameters (H0=70, Ω_m=0.3, Ω_Λ=0.7), the calculated comoving distance is {comoving_distance_gpc:.2f} Gpc.\n"
            f"3. This calculated distance is closest to option {closest_option_label} ({options[closest_option_label]} Gpc), not option {llm_final_answer_option} ({options[llm_final_answer_option]} Gpc)."
        )
        return reason

# Run the check and print the result
result = check_quasar_distance()
print(result)