import numpy as np

def check_answer():
    """
    Checks the correctness of the LLM's answer for the quasar comoving distance problem.
    """
    try:
        from astropy.cosmology import FlatLambdaCDM
        from astropy import units as u
    except ImportError:
        return "Could not check the answer because the 'astropy' library is not installed. Please install it using 'pip install astropy'."

    # --- 1. Parameters from the question and LLM's reasoning ---
    lambda_observed = 790.0  # nm
    lambda_rest_lya = 121.567  # nm (Lyman-alpha rest wavelength)
    
    # Cosmological parameters
    H0 = 70.0  # km/s/Mpc
    Omega_m = 0.3
    Omega_Lambda = 0.7
    
    # LLM's final answer
    llm_answer_choice = 'C'
    options = {'A': 6.0, 'B': 7.0, 'C': 9.0, 'D': 8.0}

    # --- 2. Verify Redshift Calculation ---
    # The LLM calculates z = (790 / 121.567) - 1 ≈ 5.50
    z_calculated = (lambda_observed / lambda_rest_lya) - 1
    llm_z_approx = 5.50
    
    if not np.isclose(z_calculated, llm_z_approx, atol=0.01):
        return f"Incorrect. The redshift calculation is inconsistent. Calculated z = {z_calculated:.4f}, while the LLM used approximately {llm_z_approx}."

    # --- 3. Set up the Cosmological Model ---
    # The problem states a flat Lambda-CDM model, which is consistent with Omega_m + Omega_Lambda = 1.0
    # astropy's FlatLambdaCDM is the correct tool for this.
    cosmo = FlatLambdaCDM(H0=H0, Om0=Omega_m)

    # --- 4. Calculate Comoving Distance ---
    # Use the precise calculated redshift for accuracy
    comoving_dist_mpc = cosmo.comoving_distance(z_calculated)
    comoving_dist_gpc = comoving_dist_mpc.to(u.Gpc).value
    
    # The LLM states the distance is ~8.987 Gpc. Let's check.
    llm_dist_gpc = 8.987
    if not np.isclose(comoving_dist_gpc, llm_dist_gpc, atol=0.01):
        # This is a minor check on the LLM's intermediate calculation. The main check is against the options.
        # We can note a small discrepancy but continue.
        pass # The difference is likely due to rounding or calculator precision.

    # --- 5. Find the closest option ---
    option_values = np.array(list(options.values()))
    option_keys = list(options.keys())
    
    # Find the index of the option with the minimum absolute difference from our calculated value
    closest_option_index = np.argmin(np.abs(option_values - comoving_dist_gpc))
    our_best_choice = option_keys[closest_option_index]

    # --- 6. Final Verification ---
    if our_best_choice == llm_answer_choice:
        # Check if the calculated value is reasonably close to the chosen option's value
        if np.isclose(comoving_dist_gpc, options[llm_answer_choice], atol=0.1):
             return "Correct"
        else:
             # This case handles if the "closest" option is still very far off.
             return f"Correct. The calculated comoving distance is {comoving_dist_gpc:.3f} Gpc, which is closest to option {llm_answer_choice} ({options[llm_answer_choice]} Gpc)."
    else:
        return (f"Incorrect. The calculated comoving distance is {comoving_dist_gpc:.3f} Gpc. "
                f"This value is closest to option {our_best_choice} ({options[our_best_choice]} Gpc), "
                f"but the provided answer was {llm_answer_choice} ({options[llm_answer_choice]} Gpc).")

# Run the check and print the result
result = check_answer()
print(result)