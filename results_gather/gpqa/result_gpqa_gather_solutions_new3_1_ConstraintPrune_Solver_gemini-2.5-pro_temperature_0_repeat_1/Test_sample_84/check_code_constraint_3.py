import math

def check_answer():
    """
    This function checks the correctness of the provided answer to the exoplanet temperature ratio problem.
    """
    # --- Given values from the problem ---
    # Planet 1
    m_p1 = 7.0  # Earth masses
    delta_lambda1 = 0.03  # Angstroms

    # Planet 2
    m_p2 = 5.0  # Earth masses
    delta_lambda2 = 0.04  # Angstroms

    # The problem states that the star's properties (T_eff, R_star) and the planets' albedos (A) are the same.
    # The specific values for the star's mass, radius, temperature, and the planets' radii are not needed for the ratio calculation.

    # --- Physics Formulas ---
    # 1. Equilibrium Temperature (T_eq)
    # T_eq is proportional to 1 / sqrt(a), where 'a' is the semi-major axis.
    # T_eq1 / T_eq2 = sqrt(a2 / a1)

    # 2. Radial Velocity Semi-amplitude (K)
    # For a circular orbit (e=0) and M_p << M_star, and assuming inclination i=90 degrees (implied by transit detection):
    # K is proportional to M_p / sqrt(a)
    # K is also directly proportional to the Doppler shift Δλ.
    # So, K1 / K2 = Δλ1 / Δλ2

    # --- Derivation ---
    # From K ∝ M_p / sqrt(a), we can write sqrt(a) ∝ M_p / K.
    # So, sqrt(a2) / sqrt(a1) = (M_p2 / K2) / (M_p1 / K1)
    # sqrt(a2 / a1) = (M_p2 / M_p1) * (K1 / K2)
    # Since T_eq1 / T_eq2 = sqrt(a2 / a1), we have:
    # T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2)

    # --- Calculation ---
    # Ratio of planet masses (Planet 2 to Planet 1)
    mass_ratio = m_p2 / m_p1

    # Ratio of RV semi-amplitudes (Planet 1 to Planet 2)
    rv_ratio = delta_lambda1 / delta_lambda2

    # Calculate the temperature ratio
    calculated_temp_ratio = mass_ratio * rv_ratio

    # --- Verification ---
    # The provided answer is 'A', which corresponds to a value of ~1.30
    llm_answer_option = 'A'
    options = {'A': 1.30, 'B': 1.05, 'C': 0.98, 'D': 0.53}
    llm_answer_value = options[llm_answer_option]

    # Check if the calculated value matches the LLM's answer value
    # We use a tolerance to account for rounding in the options.
    tolerance = 0.01
    if abs(calculated_temp_ratio - llm_answer_value) < tolerance:
        return "Correct"
    else:
        # Find the correct option
        correct_option = ''
        min_diff = float('inf')
        for option, value in options.items():
            diff = abs(calculated_temp_ratio - value)
            if diff < min_diff:
                min_diff = diff
                correct_option = option
        
        return (f"Incorrect. The provided answer is {llm_answer_option} ({llm_answer_value}), but the calculated ratio is {calculated_temp_ratio:.4f}. "
                f"The calculation is as follows:\n"
                f"1. The ratio of equilibrium temperatures T_eq1 / T_eq2 simplifies to sqrt(a2 / a1).\n"
                f"2. The ratio of semi-major axes a2 / a1 can be derived from the radial velocity formula and is equal to (M_p2 / M_p1)^2 * (K1 / K2)^2.\n"
                f"3. Combining these gives T_eq1 / T_eq2 = (M_p2 / M_p1) * (K1 / K2).\n"
                f"4. Plugging in the values: T_eq1 / T_eq2 = (5 / 7) * (0.03 / 0.04) = (5/7) * (3/4) = 15 / 28 ≈ 0.5357.\n"
                f"This value corresponds to option D (~0.53), not A (~1.30).")

# Run the check
result = check_answer()
print(result)
