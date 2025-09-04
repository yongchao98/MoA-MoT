import math

def check_correctness():
    """
    Checks the correctness of the given answer for the mean decay distance problem.
    """
    # --- Given values from the question ---
    E_X = 8.0  # Total energy in GeV
    m_X = 1.2  # Rest mass in GeV/c^2
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Constants ---
    # h-bar * c in GeV*m (from CODATA 2018 / PDG)
    HBAR_C_GEV_M = 1.97327e-16

    # --- LLM's Answer ---
    # The provided answer is B, which corresponds to the value 4.0655 * 10^-15 m.
    llm_answer_value = 4.0655e-15

    # --- Calculation ---
    # 1. Convert decay width from MeV to GeV for consistent units.
    Gamma_X_GeV = Gamma_X_MeV / 1000.0

    # 2. Check for unphysical input.
    if E_X < m_X:
        return "Incorrect. The total energy (8 GeV) cannot be less than the rest mass (1.2 GeV)."

    # 3. Calculate the particle's momentum (pc) in GeV.
    # E^2 = (pc)^2 + (mc^2)^2  =>  pc = sqrt(E^2 - (mc^2)^2)
    pc = math.sqrt(E_X**2 - m_X**2)

    # 4. Calculate the mean decay distance L in meters.
    # L = (pc / m) * (hbar*c / Gamma)
    # Note: In this formula, 'm' is the rest mass energy (m*c^2).
    calculated_distance = (pc / m_X) * (HBAR_C_GEV_M / Gamma_X_GeV)

    # --- Verification ---
    # Compare the calculated distance with the value from the LLM's answer.
    # A small relative tolerance (e.g., 0.1%) is used to account for potential
    # rounding differences in the problem's options or constants.
    if math.isclose(calculated_distance, llm_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_distance:.4e} m. "
                f"The provided answer B corresponds to a value of {llm_answer_value:.4e} m. "
                f"The calculated value does not match the answer's value within a reasonable tolerance.")

# Run the check and print the result.
result = check_correctness()
print(result)