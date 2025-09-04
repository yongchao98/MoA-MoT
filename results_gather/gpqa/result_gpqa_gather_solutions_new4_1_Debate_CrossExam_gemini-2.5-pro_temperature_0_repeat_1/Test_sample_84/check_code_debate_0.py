import math

def check_answer():
    """
    Checks the correctness of the provided answer for the exoplanet temperature ratio problem.
    """
    # Given values from the problem
    m_p1 = 7.0  # Mass of Planet 1 in Earth masses
    m_p2 = 5.0  # Mass of Planet 2 in Earth masses
    delta_lambda1 = 0.03  # Doppler shift for Planet 1 in Angstroms
    delta_lambda2 = 0.04  # Doppler shift for Planet 2 in Angstroms

    # The provided options
    options = {
        "A": 1.30,
        "B": 0.53,
        "C": 1.05,
        "D": 0.98
    }
    
    # The final answer given by the LLM
    llm_answer_key = "A"

    # --- Step-by-step calculation ---

    # The ratio of equilibrium temperatures T_eq1 / T_eq2 is proportional to sqrt(a2 / a1).
    # The semi-major axis 'a' is proportional to (M_p / K)^2, where K is the RV semi-amplitude.
    # The RV semi-amplitude K is proportional to the Doppler shift Δλ.
    # Therefore, a ∝ (M_p / Δλ)^2.
    # So, a2 / a1 = (M_p2 / Δλ2)^2 / (M_p1 / Δλ1)^2 = (M_p2 / M_p1)^2 * (Δλ1 / Δλ2)^2.
    # And T_eq1 / T_eq2 = sqrt(a2 / a1) = (M_p2 / M_p1) * (Δλ1 / Δλ2).

    # Calculate the ratio of masses (Planet 2 / Planet 1)
    mass_ratio = m_p2 / m_p1

    # Calculate the ratio of Doppler shifts (Planet 1 / Planet 2)
    doppler_shift_ratio = delta_lambda1 / delta_lambda2

    # Calculate the final temperature ratio
    calculated_ratio = mass_ratio * doppler_shift_ratio

    # Find the closest option to the calculated result
    closest_option = min(options.keys(), key=lambda k: abs(options[k] - calculated_ratio))

    # Check if the LLM's answer matches the correct option
    if llm_answer_key == closest_option:
        return "Correct"
    else:
        return (f"Incorrect. The provided answer is {llm_answer_key} ({options[llm_answer_key]}), but the calculated ratio is {calculated_ratio:.4f}. "
                f"This value is closest to option {closest_option} (~{options[closest_option]}). "
                f"The derivation is: T_eq1/T_eq2 = (M_p2/M_p1) * (Δλ₁/Δλ₂) = (5/7) * (0.03/0.04) = 15/28 ≈ 0.5357.")

# Run the check
result = check_answer()
print(result)