import math

def check_answer():
    """
    This function verifies the calculation for the exoplanet temperature ratio problem.
    """
    # Given values from the problem
    m_p1 = 7.0  # Mass of Planet 1 in Earth masses
    m_p2 = 5.0  # Mass of Planet 2 in Earth masses
    delta_lambda1 = 0.03  # Doppler shift for Planet 1 in Angstroms
    delta_lambda2 = 0.04  # Doppler shift for Planet 2 in Angstroms

    # The provided answer options
    options = {'A': 0.53, 'B': 1.05, 'C': 1.30, 'D': 0.98}
    llm_answer_key = 'A'
    llm_answer_value = options[llm_answer_key]

    # --- Step-by-step calculation ---

    # 1. The ratio of equilibrium temperatures T_eq1 / T_eq2 is sqrt(a2 / a1).
    # 2. The semi-major axis 'a' is proportional to (M_p / K)^2, where K is the RV semi-amplitude.
    # 3. The ratio of semi-major axes a2 / a1 is therefore ((M_p2 / K2) / (M_p1 / K1))^2
    #    which simplifies to (M_p2 / M_p1)^2 * (K1 / K2)^2.
    # 4. The RV semi-amplitude K is proportional to the Doppler shift Δλ.
    #    So, K1 / K2 = Δλ₁ / Δλ₂.
    # 5. Substituting this into the temperature ratio formula:
    #    T_eq1 / T_eq2 = sqrt( (M_p2 / M_p1)^2 * (Δλ₁ / Δλ₂)^2 )
    #    T_eq1 / T_eq2 = (M_p2 / M_p1) * (Δλ₁ / Δλ₂)

    # Calculate the ratio of masses (M_p2 / M_p1)
    mass_ratio = m_p2 / m_p1

    # Calculate the ratio of Doppler shifts (Δλ₁ / Δλ₂)
    doppler_shift_ratio = delta_lambda1 / delta_lambda2

    # Calculate the final temperature ratio
    temp_ratio = mass_ratio * doppler_shift_ratio

    # Check if the calculated result matches the selected option A
    # We use math.isclose to handle potential floating-point inaccuracies.
    # A tolerance of 0.01 is sufficient for this problem.
    if math.isclose(temp_ratio, llm_answer_value, rel_tol=0.01):
        return "Correct"
    else:
        # Find the correct option based on the calculation
        closest_option = min(options, key=lambda k: abs(options[k] - temp_ratio))
        
        return (f"Incorrect. The calculated ratio is {temp_ratio:.4f}. "
                f"This value is closest to option {closest_option} ({options[closest_option]}), "
                f"but the provided answer is {llm_answer_key} ({llm_answer_value}).")

# Run the check
result = check_answer()
print(result)
