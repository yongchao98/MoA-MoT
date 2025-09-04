import math

def check_correctness():
    """
    Calculates the mean decay distance of the X meson resonance and checks
    if it matches the provided answer.
    """
    # --- Given values from the question ---
    E_X_GeV = 8.0       # Production energy in GeV
    m_X_c2_GeV = 1.2    # Rest mass energy in GeV
    Gamma_X_MeV = 320.0 # Decay width in MeV

    # --- Physical constant ---
    # hbar*c is approximately 197.327 MeV*fm
    hbar_c_MeV_fm = 197.327

    # --- Convert units for consistency (use MeV) ---
    E_X_MeV = E_X_GeV * 1000
    m_X_c2_MeV = m_X_c2_GeV * 1000

    # --- Calculation Steps ---
    # 1. Calculate the momentum term, pc = sqrt(E^2 - (mc^2)^2)
    try:
        # The term inside the square root
        pc_squared_MeV2 = E_X_MeV**2 - m_X_c2_MeV**2
        if pc_squared_MeV2 < 0:
            return "Incorrect. The rest mass energy cannot be greater than the total energy."
        pc_MeV = math.sqrt(pc_squared_MeV2)
    except ValueError:
        return "Incorrect. Mathematical error during momentum calculation."

    # 2. Calculate the relativistic factor, gamma*beta = pc / mc^2
    gamma_beta = pc_MeV / m_X_c2_MeV

    # 3. Calculate the mean decay distance in femtometers (fm)
    # L = (gamma*beta) * (hbar*c / Gamma)
    L_fm = gamma_beta * (hbar_c_MeV_fm / Gamma_X_MeV)

    # 4. Convert the final result to meters (1 fm = 1e-15 m)
    calculated_L_meters = L_fm * 1e-15

    # --- Verification ---
    # The final answer from the LLM analysis is <<<B>>>.
    # The options provided in the question are:
    # A) 4.0655 * 10^-16 m
    # B) 4.0655 * 10^-15 m
    # C) 5.0223 * 10^-16 m
    # D) 5.0223 * 10^-15 m
    # Answer B corresponds to the value 4.0655 * 10^-15 m.
    expected_answer_value = 4.0655e-15

    # Check if the calculated value is close to the expected value,
    # allowing for small rounding differences in constants (e.g., 0.1% tolerance).
    if math.isclose(calculated_L_meters, expected_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_L_meters:.5e} m, "
                f"which does not match the expected answer of {expected_answer_value:.5e} m from option B.")

# Run the check
result = check_correctness()
print(result)