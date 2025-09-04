import math

def check_correctness():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the provided answer.
    """
    # --- 1. Define given values and constants ---
    # Given values from the question
    E_X_GeV = 8.0      # Total energy in GeV
    m_X_GeV = 1.2      # Rest mass energy in GeV
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # Physical constants
    # hbar * c is a useful constant in particle physics, approx. 197.327 MeV*fm
    hbar_c_MeV_fm = 197.327
    # Conversion factor from femtometers (fm) to meters (m)
    fm_to_m = 1e-15

    # --- 2. Perform the calculation ---
    # The formula for the mean decay distance (L) in the lab frame is:
    # L = (gamma * beta) * (c * tau)
    # This can be expressed in terms of given quantities as:
    # L = (sqrt(E^2 - (mc^2)^2) / mc^2) * (hbar*c / Gamma)

    # Ensure consistent units for energy (convert all to GeV)
    Gamma_X_GeV = Gamma_X_MeV / 1000.0
    hbar_c_GeV_fm = hbar_c_MeV_fm / 1000.0

    # Calculate the momentum term: pc = sqrt(E^2 - (mc^2)^2)
    # The term inside the square root must be non-negative
    energy_momentum_term = E_X_GeV**2 - m_X_GeV**2
    if energy_momentum_term < 0:
        return "Incorrect: The total energy E_X (8 GeV) must be greater than or equal to the rest mass m_X (1.2 GeV) for the particle to exist."

    pc_GeV = math.sqrt(energy_momentum_term)

    # Calculate the mean decay distance in meters
    # L = (pc / mc^2) * (hbar*c / Gamma)
    calculated_distance_m = (pc_GeV / m_X_GeV) * (hbar_c_GeV_fm / Gamma_X_GeV) * fm_to_m

    # --- 3. Compare the calculation with the provided answer ---
    # The options as listed in the question
    options = {
        'A': 5.0223e-16,
        'B': 4.0655e-16,
        'C': 5.0223e-15,
        'D': 4.0655e-15
    }

    # The final answer from the LLM to be checked
    llm_answer_letter = 'D'
    llm_answer_value = options.get(llm_answer_letter)

    if llm_answer_value is None:
        return f"Incorrect: The provided answer '{llm_answer_letter}' is not a valid option letter (A, B, C, D)."

    # Check if the calculated value is close to the value of the chosen option.
    # A relative tolerance of 0.1% is used to account for rounding of constants.
    if math.isclose(calculated_distance_m, llm_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        # If incorrect, find which option the calculation actually matches
        correct_letter = "None"
        for letter, value in options.items():
            if math.isclose(calculated_distance_m, value, rel_tol=1e-3):
                correct_letter = letter
                break
        
        reason = (f"Incorrect. The calculated mean decay distance is approximately {calculated_distance_m:.4e} m. "
                  f"This value corresponds to option '{correct_letter}' ({options.get(correct_letter):.4e} m). "
                  f"The provided answer was '{llm_answer_letter}', which corresponds to a value of {llm_answer_value:.4e} m.")
        return reason

# Run the check and print the result
print(check_correctness())