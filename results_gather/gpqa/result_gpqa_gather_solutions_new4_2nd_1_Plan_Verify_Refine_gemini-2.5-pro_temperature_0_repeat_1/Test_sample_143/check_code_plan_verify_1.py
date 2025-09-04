import math

def check_decay_distance():
    """
    Checks the correctness of the calculated mean decay distance for a meson resonance.

    The problem provides:
    - Production energy E_X = 8 GeV
    - Mass m_X = 1.2 GeV (interpreted as rest mass energy m_X * c^2)
    - Width Gamma_X = 320 MeV

    The formula for the mean decay distance L in the lab frame is:
    L = (gamma * beta) * (hbar * c / Gamma_X)
    where gamma * beta = pc / (m_X * c^2) = sqrt(E_X^2 - (m_X * c^2)^2) / (m_X * c^2)

    The final answer to check is B, which corresponds to 4.0655 * 10^-15 m.
    """

    # --- 1. Define given values and constants ---
    try:
        # Given values from the question
        E_X_GeV = 8.0
        m_X_c2_GeV = 1.2
        Gamma_X_MeV = 320.0

        # The answer to check is option B
        # The prompt lists the options as:
        # A) 5.0223 * 10^-15 m
        # B) 4.0655 * 10^-15 m
        # C) 4.0655 * 10^-16 m
        # D) 5.0223 * 10^-16 m
        # The final answer provided is <<<B>>>, so the target value is 4.0655 * 10^-15 m.
        target_answer_value = 4.0655e-15

        # Physical constant hbar * c. Using a standard precise value.
        # The minor difference between the calculated value and the option value is
        # often due to the specific precision of this constant used to generate the options.
        hbar_c_MeV_fm = 197.327  # in MeV * fm, where 1 fm = 1e-15 m

        # --- 2. Perform calculations with consistent units (MeV) ---
        # Convert GeV to MeV
        E_X_MeV = E_X_GeV * 1000
        m_X_c2_MeV = m_X_c2_GeV * 1000

        # Check if energy is sufficient for the particle to exist
        if E_X_MeV < m_X_c2_MeV:
            return "Incorrect: The total energy (8 GeV) cannot be less than the rest mass energy (1.2 GeV)."

        # Calculate the relativistic factor gamma * beta
        # pc = sqrt(E^2 - (mc^2)^2)
        pc_squared_MeV2 = E_X_MeV**2 - m_X_c2_MeV**2
        pc_MeV = math.sqrt(pc_squared_MeV2)
        
        # gamma * beta = pc / (mc^2)
        gamma_beta = pc_MeV / m_X_c2_MeV

        # Calculate the proper decay length (c * tau) in femtometers (fm)
        # c * tau = hbar * c / Gamma
        c_tau_fm = hbar_c_MeV_fm / Gamma_X_MeV

        # Calculate the mean decay distance in the lab frame (in fm)
        # L = (gamma * beta) * (c * tau)
        L_fm = gamma_beta * c_tau_fm

        # Convert the final result from fm to meters
        fm_to_m = 1e-15
        calculated_L_m = L_fm * fm_to_m

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- 3. Compare the calculated result with the target answer ---
    # We use math.isclose to account for potential rounding differences from the
    # precision of the hbar*c constant used to create the options. A relative
    # tolerance of 0.1% (1e-3) is more than enough for this.
    if math.isclose(calculated_L_m, target_answer_value, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect: The calculated mean decay distance is {calculated_L_m:.4e} m, "
                f"which does not match the expected answer of {target_answer_value:.4e} m. "
                f"The discrepancy is larger than the allowed tolerance for rounding of constants. "
                f"The calculation steps were: pc = {pc_MeV:.2f} MeV, gamma*beta = {gamma_beta:.4f}, "
                f"c*tau = {c_tau_fm:.4f} fm, leading to L = {L_fm:.4f} fm.")

# Execute the check
print(check_decay_distance())