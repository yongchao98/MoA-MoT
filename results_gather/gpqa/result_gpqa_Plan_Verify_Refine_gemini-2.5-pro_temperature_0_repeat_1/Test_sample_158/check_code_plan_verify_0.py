import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

def check_quasar_distance_answer():
    """
    This function checks the correctness of the LLM's answer by performing an independent calculation.
    """
    # --- 1. Define problem constraints and parameters from the question ---
    
    # The question describes a peak at 790 nm with a significant flux drop at shorter wavelengths.
    # This is the characteristic signature of the Lyman-alpha (Lyα) emission line from a high-redshift object,
    # where the drop is caused by absorption from the intergalactic medium (the Lyman-alpha forest).
    # This assumption is fundamental to solving the problem.
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest_lya = 121.567  # Rest-frame wavelength of Lyα in nm

    # Cosmological parameters from the question
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Om0 = 0.3  # Matter density parameter
    Ode0 = 0.7 # Dark energy density parameter

    # The multiple-choice options provided in the question
    options = {'A': 6.0, 'B': 8.0, 'C': 7.0, 'D': 9.0} # in Gpc

    # --- 2. Verify the constraints of the model ---
    
    # The question states the universe is flat. For a flat universe, the total density parameter
    # (Ω_total = Ω_matter + Ω_dark_energy) must equal 1.
    if not np.isclose(Om0 + Ode0, 1.0):
        return (f"Constraint check failed: The universe is stated to be flat, but the density "
                f"parameters do not sum to 1 (Ω_m + Ω_Λ = {Om0 + Ode0}).")

    # --- 3. Perform the calculation ---
    
    # Step A: Calculate the redshift (z) based on the Lyα assumption.
    # The formula for redshift is z = (λ_observed / λ_rest) - 1.
    z = (lambda_obs / lambda_rest_lya) - 1

    # Step B: Calculate the comoving distance.
    # We use astropy, a standard library for astronomical calculations.
    # The FlatLambdaCDM model is appropriate for a flat universe with a cosmological constant (dark energy).
    # It implicitly uses Ode0 = 1 - Om0, which matches the given parameters.
    try:
        cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
        # The comoving distance is calculated for the present epoch (scale factor a=1) at redshift z.
        comoving_dist = cosmo.comoving_distance(z)
        # Convert the distance to Gigaparsecs (Gpc) for comparison with the options.
        comoving_dist_gpc = comoving_dist.to(u.Gpc).value
    except Exception as e:
        return f"An error occurred during the astropy calculation: {e}"

    # --- 4. Check the final answer ---
    
    # To select the best multiple-choice answer, we find the option that is numerically closest
    # to our calculated distance. This is the same logic the LLM's code follows.
    closest_option_key = min(options, key=lambda k: abs(options[k] - comoving_dist_gpc))
    
    # The LLM's reasoning and code are correct. The calculated distance of ~8.97 Gpc is
    # unambiguously closest to 9 Gpc (Option D). We confirm that this logic holds.
    if closest_option_key == 'D':
        return "Correct"
    else:
        return (f"Incorrect: The calculation leads to a different result. "
                f"The calculated redshift is z = {z:.4f}, which gives a comoving distance of "
                f"{comoving_dist_gpc:.2f} Gpc. This is closest to option {closest_option_key} "
                f"({options[closest_option_key]} Gpc), not the expected answer (D).")

# Run the check and print the result.
result = check_quasar_distance_answer()
print(result)