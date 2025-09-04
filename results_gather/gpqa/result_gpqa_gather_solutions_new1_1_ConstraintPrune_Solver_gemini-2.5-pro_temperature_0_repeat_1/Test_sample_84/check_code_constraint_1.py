import math

def check_exoplanet_temperature_ratio():
    """
    This function verifies the calculation for the ratio of equilibrium temperatures
    of two exoplanets based on their mass and the radial velocity shifts they induce.
    """
    # --- Problem Data ---
    # Planet 1
    m1 = 7.0  # Mass in Earth masses
    delta_lambda1 = 0.03  # Doppler shift in Angstroms

    # Planet 2
    m2 = 5.0  # Mass in Earth masses
    delta_lambda2 = 0.04  # Doppler shift in Angstroms

    # Note: Star's mass, radius, temperature, and planets' radii are not needed
    # for the ratio calculation as they cancel out. The same albedo is a key condition.

    # --- Derivation Check ---
    # The final derived formula is: T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)
    # This is because T_eq ∝ 1/sqrt(a) and K ∝ M_p/sqrt(a), where K ∝ Δλ.
    # So, T_eq1/T_eq2 = sqrt(a2/a1) = sqrt((M_p2/K2)^2 / (M_p1/K1)^2) = (M_p2/M_p1) * (K1/K2)
    # And since K1/K2 = Δλ1/Δλ2, the formula is correct.

    # --- Calculation ---
    mass_ratio_p2_p1 = m2 / m1
    doppler_ratio_1_2 = delta_lambda1 / delta_lambda2
    
    temperature_ratio = mass_ratio_p2_p1 * doppler_ratio_1_2
    
    # --- LLM Answer Verification ---
    # The LLM's calculated value is ~0.5357 and the chosen option is D (~0.53)
    llm_calculated_value = 0.5357
    llm_option = 'D'
    
    # Check if our calculation matches the LLM's calculation
    if not math.isclose(temperature_ratio, llm_calculated_value, rel_tol=1e-3):
        return f"Calculation Mismatch: The code calculates {temperature_ratio:.4f}, but the LLM answer calculated {llm_calculated_value:.4f}."

    # Check if the LLM's chosen option is the closest one
    options = {'A': 1.05, 'B': 0.98, 'C': 1.30, 'D': 0.53}
    
    # Find the option key with the minimum absolute difference from the calculated value
    closest_option = min(options, key=lambda k: abs(options[k] - temperature_ratio))

    if closest_option != llm_option:
        return f"Incorrect Option Selection: The calculated ratio is {temperature_ratio:.4f}. The closest option is {closest_option} ({options[closest_option]}), but the LLM selected option {llm_option} ({options[llm_option]})."

    return "Correct"

# Execute the check
result = check_exoplanet_temperature_ratio()
print(result)