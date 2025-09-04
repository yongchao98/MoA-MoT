import math

def check_correctness():
    """
    This function checks the correctness of the provided answer for the mean decay distance problem.
    It recalculates the value based on fundamental physics principles and compares it to the given options.
    """
    # --- Given Parameters ---
    E_X_GeV = 8.0      # Production energy in GeV
    m_X_GeV = 1.2      # Mass in GeV/c^2
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- LLM's Answer ---
    # The LLM selected option A.
    llm_answer_key = "A"
    options = {
        "A": 4.0655e-15,
        "B": 5.0223e-16,
        "C": 5.0223e-15,
        "D": 4.0655e-16,
    }
    expected_value = options[llm_answer_key]

    # --- Physical Constants ---
    # Reduced Planck constant * speed of light (hbar*c)
    hbar_c_MeV_fm = 197.327  # Unit: MeV*fm
    # Conversion from femtometers to meters
    fm_to_m = 1e-15

    # --- Calculation ---
    # 1. Convert all energy units to MeV for consistency.
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000

    # 2. Check if the energy is sufficient for the particle to exist.
    if E_X_MeV < m_X_MeV:
        return f"Incorrect: The production energy E_X ({E_X_GeV} GeV) is less than the rest mass m_X ({m_X_GeV} GeV), which is physically impossible."

    # 3. Calculate the momentum term (pc) using the relativistic energy-momentum relation:
    # E^2 = (pc)^2 + (m*c^2)^2
    pc_squared = E_X_MeV**2 - m_X_MeV**2
    pc_MeV = math.sqrt(pc_squared)

    # 4. Calculate the mean decay distance 'd' using the formula:
    # d = (pc / m*c^2) * (hbar*c / Gamma)
    # This formula combines the Lorentz factor, velocity, and proper lifetime.
    decay_distance_fm = (pc_MeV / m_X_MeV) * (hbar_c_MeV_fm / Gamma_X_MeV)

    # 5. Convert the result from femtometers to meters.
    calculated_distance_m = decay_distance_fm * fm_to_m

    # --- Verification ---
    # Compare the calculated value with the expected value from the LLM's answer.
    # A relative tolerance of 0.1% is used to account for rounding of constants.
    if math.isclose(calculated_distance_m, expected_value, rel_tol=1e-3):
        return "Correct"
    else:
        # If the answer is incorrect, provide the reason.
        return (f"Incorrect: The calculated mean decay distance is {calculated_distance_m:.4e} m. "
                f"The LLM's answer was option {llm_answer_key}, which corresponds to a value of {expected_value:.4e} m. "
                f"The calculated value does not match the LLM's answer.")

# Run the check
result = check_correctness()
print(result)