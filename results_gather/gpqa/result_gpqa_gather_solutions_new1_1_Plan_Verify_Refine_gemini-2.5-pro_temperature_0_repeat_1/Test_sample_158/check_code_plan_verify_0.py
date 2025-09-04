import numpy as np
from scipy.integrate import quad

def check_quasar_distance():
    """
    This function verifies the comoving distance calculation for a quasar based on its spectral features.
    It follows these steps:
    1. Defines the cosmological parameters and observational data from the question.
    2. Interprets the spectral feature at 790 nm as the redshifted Lyman-alpha line.
    3. Calculates the redshift (z) of the quasar.
    4. Numerically integrates the Friedmann equation to find the comoving distance.
    5. Compares the calculated distance to the provided options to verify the given answer.
    """
    # --- Define constants and parameters from the question ---
    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest_lya = 121.6  # Rest-frame Lyman-alpha wavelength in nm

    # Cosmological parameters from the Lambda-CDM model
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Omega_m = 0.3  # Matter density parameter
    Omega_L = 0.7  # Dark energy density parameter
    # The universe is flat, so Omega_k = 0, and Omega_m + Omega_L = 1.0

    # Physical constants
    c = 299792.458  # Speed of light in km/s

    # The provided answer is 'B', which corresponds to 8 Gpc.
    expected_answer_label = 'B'
    options = {'A': 7.0, 'B': 8.0, 'C': 6.0, 'D': 9.0}

    # --- Step 1: Calculate Redshift (z) ---
    # The problem describes a peak with a flux drop at shorter wavelengths,
    # which is the classic signature of the Lyman-alpha line at high redshift.
    try:
        z = (lambda_obs / lambda_rest_lya) - 1
    except Exception as e:
        return f"Failed during redshift calculation. Error: {e}"

    # --- Step 2: Calculate Comoving Distance (D_c) ---
    # Define the integrand for the comoving distance integral.
    # E(z) = sqrt(Omega_m * (1+z)^3 + Omega_k * (1+z)^2 + Omega_L)
    # For a flat universe, Omega_k = 0.
    def integrand(z_prime):
        return 1.0 / np.sqrt(Omega_m * (1 + z_prime)**3 + Omega_L)

    # Calculate the Hubble Distance in Mpc
    D_H = c / H0

    # Perform the numerical integration from z=0 to the calculated redshift
    try:
        integral_value, integral_error = quad(integrand, 0, z)
    except Exception as e:
        return f"Failed during numerical integration. Error: {e}"

    # Calculate the comoving distance in Mpc
    Dc_Mpc = D_H * integral_value

    # Convert to Gpc for comparison with options
    Dc_Gpc = Dc_Mpc / 1000.0

    # --- Step 3: Verify the result against the options ---
    # Find the option that is numerically closest to our calculated value
    closest_option_label = min(options, key=lambda k: abs(options[k] - Dc_Gpc))

    # Check if the closest option matches the expected answer from the analysis
    if closest_option_label == expected_answer_label:
        # The calculation confirms the provided answer.
        return "Correct"
    else:
        # The calculation contradicts the provided answer.
        return (f"Incorrect. The calculation based on the problem's parameters yields a comoving distance of "
                f"{Dc_Gpc:.2f} Gpc. This value is closest to option {closest_option_label} ({options[closest_option_label]} Gpc), "
                f"not the provided answer's option {expected_answer_label} ({options[expected_answer_label]} Gpc).")

# Execute the check and print the result
result = check_quasar_distance()
print(result)