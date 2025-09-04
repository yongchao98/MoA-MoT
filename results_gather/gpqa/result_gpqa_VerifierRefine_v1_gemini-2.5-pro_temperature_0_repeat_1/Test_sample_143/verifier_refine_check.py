import math

def check_decay_distance_correctness():
    """
    Checks the correctness of the calculated mean decay distance for a meson resonance.

    The function recalculates the mean decay distance based on the given physical
    quantities and compares it to the value from the selected option.
    """
    # --- Given values from the question ---
    E_X = 8.0  # Total energy in GeV
    m_X = 1.2  # Rest mass in GeV/c^2
    Gamma_X = 320.0  # Decay width in MeV

    # --- Physical constants ---
    # The product hbar*c is a convenient constant in particle physics.
    # hbar*c ≈ 197.327 MeV·fm (where 1 fm = 10^-15 m)
    hbar_c_MeV_fm = 197.3269804

    # --- Convert units for consistency (all energies to MeV) ---
    E_X_MeV = E_X * 1000
    m_X_MeV = m_X * 1000
    Gamma_X_MeV = Gamma_X

    # --- Step 1: Calculate the Lorentz factor (gamma) ---
    # gamma = E_total / E_rest
    # In units where c=1, this is simply the ratio of energies.
    try:
        if m_X_MeV <= 0:
            return "Constraint violated: Mass must be positive."
        gamma = E_X_MeV / m_X_MeV
    except ZeroDivisionError:
        return "Error: Division by zero. Mass cannot be zero."

    # Sanity check: gamma must be >= 1 for a particle with mass
    if gamma < 1:
        return f"Constraint violated: Calculated Lorentz factor gamma ({gamma:.4f}) is less than 1, which is physically impossible."

    # --- Step 2: Calculate the velocity factor (beta) ---
    # beta = sqrt(1 - 1/gamma^2)
    beta = math.sqrt(1 - 1 / (gamma**2))

    # --- Step 3: Calculate the mean decay distance (L) ---
    # The formula is L = gamma * beta * (hbar * c) / Gamma
    # The result will be in femtometers (fm) because hbar_c is in MeV*fm
    # and Gamma_X is in MeV.
    L_fm = gamma * beta * hbar_c_MeV_fm / Gamma_X_MeV

    # --- Step 4: Convert the result to meters ---
    # 1 fm = 1e-15 m
    calculated_L_m = L_fm * 1e-15

    # --- Step 5: Compare with the provided answer ---
    # The provided answer is D) 4.0655 * 10^-15 m
    answer_value = 4.0655 * 10**-15

    # Use a relative tolerance for floating-point comparison, as is standard
    # for checking numerical results. A tolerance of 0.1% is reasonable.
    if math.isclose(calculated_L_m, answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_L_m:.4e} m, "
                f"which does not match the provided answer's value of {answer_value:.4e} m. "
                f"The calculation steps are: gamma = {gamma:.4f}, beta = {beta:.4f}, "
                f"L = {L_fm:.4f} fm = {calculated_L_m:.4e} m.")

# Execute the check and print the result
result = check_decay_distance_correctness()
print(result)