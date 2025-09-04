import numpy as np

def check_cosmology_answer():
    """
    Checks the correctness of the provided answer by recalculating the comoving distance.
    """
    try:
        from astropy.cosmology import FlatLambdaCDM
        import astropy.units as u
    except ImportError:
        return "The 'astropy' library is required to run this check. Please install it using 'pip install astropy'."

    # --- Step 1: Define parameters from the question ---
    
    # Observational data
    lambda_obs = 790.0  # Observed wavelength in nm
    lambda_rest_lya = 121.567  # Rest-frame wavelength of Lyman-alpha in nm
    
    # Cosmological parameters
    H0 = 70.0  # Hubble constant in km/s/Mpc
    Om0 = 0.3   # Matter density parameter
    # For a flat universe, Ode0 (dark energy density) is 1 - Om0 = 0.7
    
    # Multiple choice options from the question
    options = {
        'A': 9.0,
        'B': 6.0,
        'C': 8.0,
        'D': 7.0
    }
    
    # The final answer provided by the LLM to be checked
    llm_answer_choice = 'A'
    
    # --- Step 2: Calculate the redshift (z) ---
    z = (lambda_obs / lambda_rest_lya) - 1
    
    # --- Step 3: Calculate the comoving distance using astropy ---
    
    # Define the cosmological model
    cosmo = FlatLambdaCDM(H0=H0, Om0=Om0)
    
    # Calculate the comoving distance
    comoving_dist_mpc = cosmo.comoving_distance(z)
    
    # Convert the result to Gigaparsecs (Gpc) to match the options
    calculated_dist_gpc = comoving_dist_mpc.to(u.Gpc).value
    
    # --- Step 4: Verify the answer ---
    
    # Find which option is numerically closest to our calculated value
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_dist_gpc))
    
    # Check if the LLM's answer matches the closest option
    if llm_answer_choice == closest_option_key:
        return "Correct"
    else:
        reason = (
            f"The provided answer '{llm_answer_choice}' ({options[llm_answer_choice]} Gpc) is incorrect.\n"
            f"1. The calculated redshift is z ≈ {z:.3f}.\n"
            f"2. Using a standard cosmological calculation (astropy), the comoving distance is ≈ {calculated_dist_gpc:.3f} Gpc.\n"
            f"3. This value is closest to option '{closest_option_key}' ({options[closest_option_key]} Gpc), not '{llm_answer_choice}'."
        )
        return reason

# Run the check
result = check_cosmology_answer()
print(result)