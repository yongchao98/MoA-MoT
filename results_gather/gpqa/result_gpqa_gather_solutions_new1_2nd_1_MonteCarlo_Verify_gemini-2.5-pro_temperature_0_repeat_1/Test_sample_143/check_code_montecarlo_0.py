import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided answer for the physics problem.
    It calculates the mean decay distance from the given parameters and compares it
    to the selected option.
    """

    # --- Define problem constraints and given values ---
    E_X = 8.0  # Total energy in GeV
    m_X = 1.2  # Rest mass energy (m_X * c^2) in GeV
    Gamma_X = 320.0  # Decay width in MeV

    # Physical constant hbar*c
    # Using a standard value often used in particle physics problems.
    # The minor variations in this constant can lead to small differences in the final digits.
    hbar_c_MeV_fm = 197.327  # in MeV * fm

    # The options provided in the question
    options = {
        'A': 5.0223e-16,
        'B': 5.0223e-15,
        'C': 4.0655e-15,
        'D': 4.0655e-16
    }

    # The answer provided by the LLM to be checked
    llm_answer_key = 'C'
    llm_answer_value = options[llm_answer_key]

    # --- Perform the calculation from first principles ---

    # The formula for the mean decay distance L in the lab frame is:
    # L = (sqrt(E_X^2 - (m_X*c^2)^2) / (m_X*c^2)) * (hbar*c / Gamma_X)

    # Step 1: Ensure all energy units are consistent. We will use GeV.
    Gamma_X_GeV = Gamma_X / 1000.0  # Convert MeV to GeV

    # Step 2: Calculate the momentum term (pc) in GeV
    # pc = sqrt(E^2 - (mc^2)^2)
    try:
        pc_GeV = math.sqrt(E_X**2 - m_X**2)
    except ValueError:
        return "Calculation Error: The total energy (E_X) must be greater than the rest mass energy (m_X)."

    # Step 3: Calculate the kinematic factor (gamma * beta)
    # This is a dimensionless quantity.
    kinematic_factor = pc_GeV / m_X

    # Step 4: Calculate the decay length factor (hbar*c / Gamma_X)
    # To get the result in meters, we need to convert hbar*c and units appropriately.
    # hbar*c = 197.327 MeV*fm = 0.197327 GeV*fm
    # 1 fm = 1e-15 m
    hbar_c_GeV_m = 0.197327 * 1e-15 # in GeV * m
    
    decay_length_factor_m = hbar_c_GeV_m / Gamma_X_GeV

    # Step 5: Calculate the final mean decay distance L in meters
    calculated_L_m = kinematic_factor * decay_length_factor_m

    # --- Verify the correctness of the LLM's answer ---

    # Check if the LLM's chosen option matches the calculated value.
    # We use math.isclose() to handle potential floating-point inaccuracies.
    # A relative tolerance of 0.1% (1e-3) is reasonable for this kind of physics problem
    # due to rounding of physical constants.
    if not math.isclose(calculated_L_m, llm_answer_value, rel_tol=1e-3):
        # Find the closest option to our calculation
        min_diff = float('inf')
        best_match_key = None
        for key, value in options.items():
            diff = abs(calculated_L_m - value)
            if diff < min_diff:
                min_diff = diff
                best_match_key = key
        
        return (f"Incorrect. The calculated mean decay distance is approximately {calculated_L_m:.4e} m. "
                f"This value is closest to option {best_match_key} ({options[best_match_key]:.4e} m). "
                f"The provided answer was option {llm_answer_key} ({llm_answer_value:.4e} m), which does not match the calculation.")

    return "Correct"

# Run the check
result = check_correctness_of_answer()
print(result)