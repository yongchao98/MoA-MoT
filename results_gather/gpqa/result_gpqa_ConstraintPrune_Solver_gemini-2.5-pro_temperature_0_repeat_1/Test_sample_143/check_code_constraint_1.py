import math

def check_decay_distance():
    """
    Calculates the mean decay distance of a meson resonance and checks it against a given answer.
    """
    # --- Given values from the question ---
    # Production energy in GeV
    E_X_GeV = 8.0
    # Mass in GeV/c^2
    m_X_GeV = 1.2
    # Width in MeV
    Gamma_X_MeV = 320.0

    # --- Constants ---
    # h-bar * c in MeV*fm (1 fm = 1e-15 m)
    hbar_c_MeV_fm = 197.327

    # --- The answer to check (Option C) ---
    # Convert meters to femtometers for comparison
    provided_answer_m = 4.0655e-15
    provided_answer_fm = provided_answer_m * 1e15

    # --- Calculation ---

    # 1. Calculate the Lorentz factor (gamma)
    # gamma = E / (m*c^2). Since E and m are in energy units, it's a direct ratio.
    gamma = E_X_GeV / m_X_GeV

    # 2. Calculate the velocity factor (beta = v/c)
    # beta = sqrt(1 - 1/gamma^2)
    try:
        beta = math.sqrt(1 - 1 / gamma**2)
    except ValueError:
        return "Error: Cannot calculate beta. This implies E < m, which is physically impossible for a real particle."

    # 3. Calculate the mean decay distance (d)
    # Formula: d = beta * gamma * (hbar*c / Gamma)
    # The result will be in femtometers (fm) because hbar_c is in MeV*fm and Gamma is in MeV.
    mean_decay_distance_fm = beta * gamma * (hbar_c_MeV_fm / Gamma_X_MeV)

    # --- Verification ---

    # Check if the calculated distance is close to the provided answer's value
    # Use a relative tolerance of 0.1% for floating point comparison
    if math.isclose(mean_decay_distance_fm, provided_answer_fm, rel_tol=1e-3):
        return "Correct"
    else:
        # Convert calculated distance back to meters for the error message
        calculated_distance_m = mean_decay_distance_fm / 1e15
        return (
            f"Incorrect. The provided answer is {provided_answer_m:.4e} m, but the calculated value is {calculated_distance_m:.4e} m.\n"
            f"Constraint check:\n"
            f"- Lorentz factor (gamma) = E / m = {E_X_GeV} GeV / {m_X_GeV} GeV = {gamma:.4f}\n"
            f"- Velocity factor (beta) = sqrt(1 - 1/gamma^2) = {beta:.4f}\n"
            f"- Mean decay distance = beta * gamma * (hbar*c / Gamma) = {beta:.4f} * {gamma:.4f} * ({hbar_c_MeV_fm} MeV*fm / {Gamma_X_MeV} MeV) = {mean_decay_distance_fm:.4f} fm\n"
            f"- This corresponds to {calculated_distance_m:.4e} m, which does not match the answer from option C."
        )

# Run the check
result = check_decay_distance()
print(result)