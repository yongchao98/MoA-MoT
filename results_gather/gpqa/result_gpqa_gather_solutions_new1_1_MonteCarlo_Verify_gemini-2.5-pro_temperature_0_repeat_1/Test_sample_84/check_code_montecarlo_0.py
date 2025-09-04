import math

def check_answer():
    """
    This function checks the correctness of the answer to the exoplanet temperature ratio problem.
    """
    # --- Given values from the question ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # --- Physics Derivation ---
    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is given by:
    # T_eq1 / T_eq2 = sqrt(a2 / a1)
    # where a1 and a2 are the semi-major axes.

    # The radial velocity semi-amplitude K is proportional to M_p / sqrt(a).
    # K ∝ M_p / sqrt(a)
    # Also, K is proportional to the Doppler shift Δλ.
    # So, Δλ ∝ M_p / sqrt(a)
    # Rearranging for sqrt(a): sqrt(a) ∝ M_p / Δλ

    # Substitute this into the temperature ratio equation:
    # T_eq1 / T_eq2 = sqrt(a2) / sqrt(a1)
    # T_eq1 / T_eq2 = (M_p2 / Δλ2) / (M_p1 / Δλ1)
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ1 / Δλ2)

    # --- Calculation ---
    mass_ratio = m_p2 / m_p1
    shift_ratio = delta_lambda1 / delta_lambda2
    temperature_ratio = mass_ratio * shift_ratio

    # --- Options from the question ---
    options = {
        "A": 0.53,
        "B": 1.30,
        "C": 0.98,
        "D": 1.05
    }
    
    # The final answer provided by the LLM
    llm_final_answer_key = 'A'

    # --- Verification ---
    # Find the option that is numerically closest to the calculated result.
    closest_option = None
    min_difference = float('inf')

    for key, value in options.items():
        difference = abs(temperature_ratio - value)
        if difference < min_difference:
            min_difference = difference
            closest_option = key

    # Check if the LLM's answer matches the calculated correct option.
    if llm_final_answer_key == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The calculated ratio is (5/7) * (0.03/0.04) = 15/28 ≈ {temperature_ratio:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}). "
                f"The provided answer was {llm_final_answer_key} ({options[llm_final_answer_key]}).")

# Execute the check and print the result
result = check_answer()
print(result)