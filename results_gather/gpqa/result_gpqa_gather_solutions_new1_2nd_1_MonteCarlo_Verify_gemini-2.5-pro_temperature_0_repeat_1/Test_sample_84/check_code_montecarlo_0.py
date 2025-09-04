import math

def check_correctness():
    """
    Checks the correctness of the provided answer for the exoplanet temperature ratio problem.

    The key physical relationships are:
    1. Equilibrium Temperature (T_eq) is proportional to 1 / sqrt(a), where 'a' is the semi-major axis.
       Therefore, T_eq1 / T_eq2 = sqrt(a2 / a1).
    2. Radial Velocity semi-amplitude (K) is proportional to M_p / sqrt(a), where M_p is the planet's mass.
       Therefore, sqrt(a) is proportional to M_p / K.
    3. Combining these gives: T_eq1 / T_eq2 = (M_p2 / K2) / (M_p1 / K1) = (M_p2 / M_p1) * (K1 / K2).
    4. The Doppler shift (Δλ) is directly proportional to K, so K1 / K2 = Δλ1 / Δλ2.
    5. The final formula is: T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2).
    """
    
    # --- Problem Data ---
    # Planet 1
    M_p1 = 7  # Earth masses
    delta_lambda1 = 0.03  # Angstroms
    
    # Planet 2
    M_p2 = 5  # Earth masses
    delta_lambda2 = 0.04  # Angstroms
    
    # --- Options ---
    options = {
        "A": 0.98,
        "B": 0.53,
        "C": 1.30,
        "D": 1.05
    }
    
    # --- LLM's Answer ---
    llm_answer_choice = "B"
    
    # --- Calculation ---
    # Calculate the ratio of masses (M_p2 / M_p1)
    mass_ratio = M_p2 / M_p1
    
    # Calculate the ratio of Doppler shifts (Δλ1 / Δλ2)
    doppler_shift_ratio = delta_lambda1 / delta_lambda2
    
    # Calculate the final temperature ratio
    calculated_temp_ratio = mass_ratio * doppler_shift_ratio
    
    # --- Verification ---
    # Check if the LLM's chosen option is the closest to the calculated value
    
    # Find the option that is numerically closest to the calculated result
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_temp_ratio))
    
    if llm_answer_choice == closest_option:
        # Further check if the numerical value is consistent
        expected_value = 15 / 28
        if math.isclose(calculated_temp_ratio, expected_value, rel_tol=1e-9):
            return "Correct"
        else:
            return f"Incorrect. The logic is correct, but there's a calculation error. Expected {expected_value}, but calculated {calculated_temp_ratio}."
    else:
        return (f"Incorrect. The calculated temperature ratio is {calculated_temp_ratio:.4f}. "
                f"This value is closest to option {closest_option} (~{options[closest_option]}), "
                f"but the provided answer was option {llm_answer_choice} (~{options[llm_answer_choice]}).")

# Run the check
result = check_correctness()
print(result)