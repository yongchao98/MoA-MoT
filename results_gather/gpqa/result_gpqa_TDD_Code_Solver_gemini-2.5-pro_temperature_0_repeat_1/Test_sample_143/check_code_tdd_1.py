import math

def check_decay_distance():
    """
    This function calculates the mean decay distance of a meson resonance
    based on its energy, mass, and decay width, and verifies it against a given answer.
    """
    # --- Physical Constants ---
    # hbar * c in MeV*fm (from Particle Data Group)
    HBAR_C_MEV_FM = 197.3269804

    # --- Given values from the question ---
    E_X_GeV = 8.0      # Production energy in GeV
    m_X_GeV = 1.2      # Mass in GeV/c^2
    Gamma_X_MeV = 320.0  # Decay width in MeV

    # --- Answer to check (Option B) ---
    # The value from option B is 4.0655 * 10^-15 m
    answer_b_m = 4.0655e-15

    # --- Calculation ---

    # 1. Convert all energy units to MeV for consistency.
    E_X_MeV = E_X_GeV * 1000
    m_X_c2_MeV = m_X_GeV * 1000

    # 2. Calculate the momentum (pc) of the particle in MeV.
    # From the relativistic energy-momentum relation: E^2 = (pc)^2 + (mc^2)^2
    # So, pc = sqrt(E^2 - (mc^2)^2)
    try:
        # Ensure energy is greater than rest mass energy
        if E_X_MeV <= m_X_c2_MeV:
            return "Constraint not satisfied: The production energy E_X must be greater than the rest mass energy m_X*c^2."
        pc_MeV = math.sqrt(E_X_MeV**2 - m_X_c2_MeV**2)
    except ValueError:
        # This case is handled above, but included for robustness
        return "Calculation Error: Cannot take the square root of a negative number. E_X must be > m_X*c^2."

    # 3. Calculate the rest-frame mean decay length (c*tau).
    # c*tau = (hbar*c) / Gamma
    c_tau_fm = HBAR_C_MEV_FM / Gamma_X_MeV

    # 4. Calculate the lab-frame mean decay distance (d).
    # A direct formula is d = (pc / (m*c^2)) * (c*tau).
    d_fm = (pc_MeV / m_X_c2_MeV) * c_tau_fm

    # 5. Convert the final distance from femtometers (fm) to meters (m).
    # 1 fm = 1e-15 m
    calculated_d_m = d_fm * 1e-15

    # --- Verification ---
    # Check if the calculated value is close to the provided answer from option B.
    # A relative tolerance of 0.1% (1e-3) is appropriate to account for
    # differences in the precision of constants used.
    if math.isclose(calculated_d_m, answer_b_m, rel_tol=1e-3):
        return "Correct"
    else:
        return (f"Incorrect. The calculated mean decay distance is {calculated_d_m:.4e} m, "
                f"which does not match the provided answer of {answer_b_m:.4e} m. "
                f"The calculated value is approximately {d_fm:.4f} fm, while the answer corresponds to {answer_b_m / 1e-15:.4f} fm.")

# Run the check and print the result.
result = check_decay_distance()
print(result)