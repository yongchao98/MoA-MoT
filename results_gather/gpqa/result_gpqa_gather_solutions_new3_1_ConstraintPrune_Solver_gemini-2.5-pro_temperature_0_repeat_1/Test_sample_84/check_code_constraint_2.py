import math

def check_exoplanet_temperature_ratio():
    """
    This function verifies the calculation for the ratio of equilibrium temperatures
    of two exoplanets based on the provided data.
    """
    # --- Given Data ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- LLM's Answer ---
    llm_answer_choice = 'D'
    options = {
        'A': 1.30,
        'B': 1.05,
        'C': 0.98,
        'D': 0.53
    }

    # --- Derivation and Calculation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 simplifies to:
    # (M_p2 / M_p1) * (K1 / K2)
    # where K is the radial velocity semi-amplitude.
    # Since K is proportional to the Doppler shift Δλ, the ratio K1/K2 is equal to Δλ1/Δλ2.
    # So, T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    mass_ratio = m_p2 / m_p1
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    
    calculated_temperature_ratio = mass_ratio * doppler_shift_ratio
    
    # --- Verification ---
    # Check if the calculation is correct
    expected_value = 15.0 / 28.0
    if not math.isclose(calculated_temperature_ratio, expected_value, rel_tol=1e-5):
        return (f"Calculation is incorrect. "
                f"Expected T_eq1/T_eq2 = (5/7) * (0.03/0.04) = {expected_value:.4f}. "
                f"The code calculated {calculated_temperature_ratio:.4f}.")

    # Check if the chosen option matches the result
    llm_answer_value = options.get(llm_answer_choice)
    if llm_answer_value is None:
        return f"The answer choice '{llm_answer_choice}' is not a valid option."

    # Find the option that is numerically closest to the calculated result
    closest_option_key = min(options, key=lambda k: abs(options[k] - calculated_temperature_ratio))

    if llm_answer_choice != closest_option_key:
        return (f"The final answer choice is incorrect. "
                f"The calculated ratio is {calculated_temperature_ratio:.4f}, "
                f"which is closest to option {closest_option_key} ({options[closest_option_key]}), "
                f"not option {llm_answer_choice} ({llm_answer_value}).")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
print(check_exoplanet_temperature_ratio())