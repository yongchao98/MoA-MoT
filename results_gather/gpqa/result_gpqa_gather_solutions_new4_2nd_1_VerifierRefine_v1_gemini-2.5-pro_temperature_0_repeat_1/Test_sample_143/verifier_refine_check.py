import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against the given answer.
    """
    # --- Given values from the question ---
    E_X = 8.0  # Total energy in GeV
    m_X_c2 = 1.2  # Rest mass energy in GeV
    Gamma_X = 320.0  # Decay width in MeV

    # --- Physical constants ---
    # Using the standard value for h-bar * c in MeV*fm simplifies the calculation
    hbar_c_MeV_fm = 197.327  # in MeV * fm
    # Conversion factor from femtometers to meters
    fm_to_m = 1e-15

    # --- The answer to be checked ---
    # The final answer provided is 'A', which corresponds to 4.0655 * 10^-15 m
    # from the given list of options.
    expected_value = 4.0655e-15

    # --- Step 1: Calculate the relativistic factor (gamma * beta) ---
    # This term is dimensionless.
    # Ensure the energy is greater than the rest mass energy to avoid math domain errors.
    if E_X <= m_X_c2:
        return "Constraint not satisfied: Total energy (E_X) must be greater than rest mass energy (m_X_c2)."
    
    pc_squared = E_X**2 - m_X_c2**2
    pc = math.sqrt(pc_squared)
    gamma_beta = pc / m_X_c2

    # --- Step 2: Calculate the proper decay length (c * tau) ---
    # This is equal to (hbar * c) / Gamma.
    # The units (MeV*fm / MeV) cancel out, leaving the result in femtometers (fm).
    proper_decay_length_fm = hbar_c_MeV_fm / Gamma_X

    # --- Step 3: Calculate the mean decay distance in the lab frame ---
    # L = (gamma * beta) * (proper decay length)
    L_fm = gamma_beta * proper_decay_length_fm

    # --- Step 4: Convert the final result to meters ---
    calculated_L_m = L_fm * fm_to_m

    # --- Step 5: Compare the calculated value with the expected answer ---
    # Use math.isclose() for robust floating-point comparison, allowing for minor
    # differences due to the precision of physical constants. A relative tolerance
    # of 0.1% (1e-3) is appropriate.
    if math.isclose(calculated_L_m, expected_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is approximately {calculated_L_m:.4e} m. "
                f"The provided answer value is {expected_value:.4e} m. The values do not match.")

# Run the check
result = check_decay_distance()
print(result)