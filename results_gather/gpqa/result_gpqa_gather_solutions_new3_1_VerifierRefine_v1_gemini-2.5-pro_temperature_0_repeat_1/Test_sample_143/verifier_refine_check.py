import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the provided answer.
    """
    # Given parameters from the question
    E_X_GeV = 8.0  # Production energy in GeV
    m_X_GeV = 1.2  # Mass in GeV (implicitly GeV/c^2, so this is rest energy)
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # Physical constant
    # hbar * c ≈ 197.327 MeV*fm
    hbar_c_MeV_fm = 197.327

    # --- Calculation ---
    # For consistency, convert all energy units to MeV
    E_X_MeV = E_X_GeV * 1000
    m_X_MeV = m_X_GeV * 1000

    # The formula for the mean decay distance L is:
    # L = gamma * beta * c * tau
    # where tau = hbar / Gamma
    # So, L = (gamma * beta) * (hbar * c / Gamma)

    # We can calculate the term (gamma * beta) from the energy-momentum relation:
    # E^2 = (pc)^2 + (mc^2)^2  =>  pc = sqrt(E^2 - (mc^2)^2)
    # And since pc = gamma * beta * mc^2:
    # gamma * beta = pc / mc^2 = sqrt(E^2 - (mc^2)^2) / mc^2

    try:
        # Calculate the numerator of gamma*beta: pc
        pc_MeV = math.sqrt(E_X_MeV**2 - m_X_MeV**2)

        # Calculate gamma*beta
        gamma_beta = pc_MeV / m_X_MeV

        # Calculate the mean decay distance in femtometers (fm)
        L_fm = gamma_beta * (hbar_c_MeV_fm / Gamma_X_MeV)

        # Convert the result to meters (1 fm = 1e-15 m)
        calculated_L_m = L_fm * 1e-15
    except ValueError:
        return "Incorrect: Calculation failed. The energy must be greater than the mass for a real particle."

    # The final answer provided by the LLM is 'A', which corresponds to 4.0655 * 10^-15 m.
    # Let's check if our calculated value matches this.
    expected_answer_value = 4.0655e-15

    # We use math.isclose() to handle potential floating-point inaccuracies.
    # A relative tolerance of 0.1% (1e-3) is reasonable for this kind of physics problem
    # due to slight variations in constants used.
    if math.isclose(calculated_L_m, expected_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect: The calculated mean decay distance is approximately {calculated_L_m:.4e} m. "
                f"The provided answer is {expected_answer_value:.4e} m. The values do not match.")

# Run the check
result = check_decay_distance()
print(result)