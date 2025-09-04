import math

def check_correctness():
    """
    Checks the correctness of the provided answer for the mean decay distance problem.

    The mean decay distance L is calculated using the formula:
    L = (pc / mc^2) * (hbar*c / Gamma)
    where pc = sqrt(E^2 - (mc^2)^2).
    """

    # --- Given values from the question ---
    E_X_GeV = 8.0  # Total energy in GeV
    m_X_GeV = 1.2  # Rest mass energy in GeV
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Physical constants ---
    # hbar * c in MeV*fm. 1 fm = 1e-15 m. Using a high-precision value.
    hbar_c_MeV_fm = 197.3269804

    # --- Unit Conversion ---
    # Convert all energy units to MeV for consistency.
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000

    # --- Calculation ---
    # 1. Calculate pc = sqrt(E^2 - (mc^2)^2)
    # This term represents the momentum multiplied by the speed of light.
    try:
        pc_squared = E_X_MeV**2 - m_X_MeV**2
        if pc_squared < 0:
            return "Calculation Error: The total energy (E) must be greater than or equal to the rest mass energy (m*c^2). E^2 - (mc^2)^2 is negative."
        pc_MeV = math.sqrt(pc_squared)
    except ValueError:
        return "Calculation Error: Could not compute the square root of E^2 - (mc^2)^2."

    # 2. Calculate the kinematic factor (pc / mc^2)
    # This is equivalent to gamma * beta.
    kinematic_factor = pc_MeV / m_X_MeV

    # 3. Calculate the decay length factor (hbar*c / Gamma)
    # This represents the proper decay length (c*tau).
    decay_factor_fm = hbar_c_MeV_fm / Gamma_X_MeV

    # 4. Calculate the mean decay distance in femtometers (fm)
    L_fm = kinematic_factor * decay_factor_fm

    # 5. Convert the final result to meters
    calculated_distance_m = L_fm * 1e-15

    # --- Verification ---
    # The provided answer from the LLM is 'C'.
    llm_answer_choice = 'C'

    # The options given in the question.
    options = {
        'A': 4.0655e-16,
        'B': 5.0223e-15,
        'C': 4.0655e-15,
        'D': 5.0223e-16
    }

    # Check if the LLM's chosen option is valid.
    if llm_answer_choice not in options:
        return f"Invalid Answer: The final answer choice '{llm_answer_choice}' is not one of the valid options (A, B, C, D)."

    # Get the value corresponding to the LLM's answer.
    llm_answer_value = options[llm_answer_choice]

    # Compare the calculated value with the LLM's answer value.
    # A small relative tolerance is used to account for potential rounding differences
    # in the physical constants used to generate the problem's options.
    tolerance = 0.001  # 0.1%
    if math.isclose(calculated_distance_m, llm_answer_value, rel_tol=tolerance):
        return "Correct"
    else:
        # If the answer is incorrect, find the correct option.
        correct_option = None
        for option, value in options.items():
            if math.isclose(calculated_distance_m, value, rel_tol=tolerance):
                correct_option = option
                break
        
        reason = (f"Incorrect. The calculated mean decay distance is approximately {calculated_distance_m:.4e} m. "
                  f"The provided answer was '{llm_answer_choice}', which corresponds to a value of {llm_answer_value:.4e} m. ")
        
        if correct_option:
            reason += f"The calculated value actually matches option '{correct_option}'."
        else:
            reason += "The calculated value does not match any of the given options within the tolerance."
            
        return reason

# Run the check and print the result.
print(check_correctness())