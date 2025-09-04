import numpy as np
from scipy.integrate import quad

def check_cosmology_answer():
    """
    Checks the correctness of the provided answer for the quasar comoving distance problem.

    The function performs the following steps:
    1. Defines the physical and cosmological constants from the question.
    2. Calculates the redshift (z) based on the observed Lyman-alpha line.
    3. Numerically integrates the Friedmann equation to find the comoving distance.
    4. Compares the calculated distance to the provided options to find the best match.
    5. Checks if the best match corresponds to the given answer.
    """
    # Step 1: Define constants and inputs from the question
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest = 121.567  # Rest-frame wavelength of Lyman-alpha in nm
    
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter
    c = 299792.458  # Speed of light in km/s

    # The options as provided in the final analysis block
    options = {'A': 6, 'B': 7, 'C': 9, 'D': 8}
    # The final answer provided by the LLM
    given_answer_letter = 'C'

    # Step 2: Calculate the redshift (z)
    # The formula is z = (lambda_obs / lambda_rest) - 1
    calculated_z = (lambda_obs / lambda_rest) - 1
    
    # The answer's reasoning states z is approx 5.5. Let's check this.
    if not np.isclose(calculated_z, 5.5, atol=0.01):
        return f"Incorrect redshift calculation in reasoning. Expected z≈5.5, but calculated z={calculated_z:.3f}."

    # Step 3: Calculate the comoving distance
    # The comoving distance D_c is given by the integral:
    # D_c = (c / H0) * integral from 0 to z of [dz' / sqrt(Omega_m*(1+z')^3 + Omega_L)]
    
    # Define the function to be integrated
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Perform the numerical integration
    try:
        integral_result, _ = quad(integrand, 0, calculated_z)
    except Exception as e:
        return f"An error occurred during numerical integration: {e}"

    # Calculate the Hubble distance in Mpc
    D_H = c / H0
    
    # Calculate the comoving distance in Mpc and then convert to Gpc
    Dc_Mpc = D_H * integral_result
    Dc_Gpc = Dc_Mpc / 1000.0

    # The answer's reasoning states the distance is ~8.99 Gpc. Let's check this.
    if not np.isclose(Dc_Gpc, 8.99, atol=0.1):
        return f"Incorrect comoving distance calculation. The reasoning states ~8.99 Gpc, but the calculated value is {Dc_Gpc:.2f} Gpc."

    # Step 4: Find the best matching option
    # Find the option key (letter) that has the minimum absolute difference with the calculated distance
    best_option_letter = min(options, key=lambda k: abs(options[k] - Dc_Gpc))
    
    # Step 5: Check if the best match corresponds to the given answer
    if best_option_letter == given_answer_letter:
        return "Correct"
    else:
        return (f"Incorrect final answer. The calculated comoving distance is {Dc_Gpc:.2f} Gpc. "
                f"This is closest to option {best_option_letter} ({options[best_option_letter]} Gpc), "
                f"but the provided answer was {given_answer_letter} ({options[given_answer_letter]} Gpc).")

# Run the check
result = check_cosmology_answer()
print(result)