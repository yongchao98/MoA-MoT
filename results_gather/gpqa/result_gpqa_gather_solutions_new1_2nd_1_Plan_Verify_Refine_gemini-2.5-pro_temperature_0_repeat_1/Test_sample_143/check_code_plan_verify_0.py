import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the provided answer.
    """
    # --- Given Physical Quantities ---
    E_X = 8.0  # Total energy in GeV
    m_X_c2 = 1.2  # Rest mass energy in GeV
    Gamma_X = 320.0  # Decay width in MeV

    # --- Physical Constant ---
    # hbar*c is approximately 197.327 MeV·fm.
    # 1 fm (femtometer) = 1e-15 m.
    hbar_c_MeV_fm = 197.327

    # --- Unit Conversions for Consistency ---
    # Convert all energy units to MeV for calculation with hbar_c.
    E_X_MeV = E_X * 1000
    m_X_c2_MeV = m_X_c2 * 1000
    # Gamma_X is already in MeV.

    # --- Calculation ---
    # The formula for the mean decay distance (L) in the lab frame is:
    # L = (sqrt(E_X^2 - (m_X*c^2)^2) / (m_X*c^2)) * (hbar*c / Gamma_X)
    # This is equivalent to L = (pc / mc^2) * (hbar*c / Gamma)

    # 1. Calculate the momentum term (pc) in MeV
    try:
        pc_squared_MeV2 = E_X_MeV**2 - m_X_c2_MeV**2
        if pc_squared_MeV2 < 0:
            return "Incorrect: The energy must be greater than the rest mass energy. Cannot calculate sqrt of a negative number."
        pc_MeV = math.sqrt(pc_squared_MeV2)
    except ValueError:
        return "Incorrect: Mathematical error during momentum calculation."

    # 2. Calculate the kinematic factor (dimensionless)
    kinematic_factor = pc_MeV / m_X_c2_MeV

    # 3. Calculate the decay length factor in femtometers (fm)
    decay_factor_fm = hbar_c_MeV_fm / Gamma_X

    # 4. Calculate the mean decay distance in femtometers (fm)
    L_fm = kinematic_factor * decay_factor_fm

    # 5. Convert the final result to meters
    L_m = L_fm * 1e-15

    # --- Verification ---
    # The provided answer is C, which corresponds to 4.0655 * 10^-15 m.
    expected_answer_value = 4.0655e-15
    
    # Check if the calculated value is close to the expected answer.
    # A relative tolerance of 0.1% is reasonable due to precision of physical constants.
    if math.isclose(L_m, expected_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is approximately {L_m:.4e} m. "
                f"The value from the selected answer C is {expected_answer_value:.4e} m. "
                f"These values do not match.")

# Run the check and print the result
result = check_decay_distance()
print(result)