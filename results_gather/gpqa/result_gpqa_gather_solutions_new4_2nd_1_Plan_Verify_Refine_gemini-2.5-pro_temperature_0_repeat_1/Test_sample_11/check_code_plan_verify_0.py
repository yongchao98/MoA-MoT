import math

def check_physics_decay_answer():
    """
    Checks the correctness of the calculated kinetic energies for the pion decay.
    Pi(+) -> mu(+) + nu
    """
    # --- Given Information from the Question ---
    m_pi_c2 = 139.6  # Rest mass energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest mass energy of mu(+) in MeV

    # --- The Answer to Verify ---
    # The provided answer is C, which corresponds to the values 4.12 MeV and 29.8 MeV.
    # We check if our calculation matches these values.
    expected_ke_muon = 4.12
    expected_ke_neutrino = 29.8

    # --- Calculation based on Physics Principles ---
    # Calculate the kinetic energy of the neutrino (KE_nu) using the derived formula.
    # For a massless particle like the neutrino, its total energy is its kinetic energy.
    try:
        ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation Error: The rest mass of the pion (m_pi_c2) cannot be zero."

    # The total kinetic energy released (Q-value) is the difference in rest mass energy.
    ke_total = m_pi_c2 - m_mu_c2

    # The kinetic energy of the muon is the remainder of the total kinetic energy.
    ke_mu = ke_total - ke_nu

    # --- Verification Step ---
    # Compare the calculated values with the expected values from the answer.
    # We use math.isclose() to handle potential floating-point inaccuracies.
    # The tolerance is set based on the precision of the options.
    is_muon_ke_correct = math.isclose(ke_mu, expected_ke_muon, abs_tol=0.01)
    is_neutrino_ke_correct = math.isclose(ke_nu, expected_ke_neutrino, abs_tol=0.01)

    if is_muon_ke_correct and is_neutrino_ke_correct:
        return "Correct"
    else:
        # If the answer is incorrect, provide the calculated values for diagnosis.
        reason = "Incorrect. The calculated kinetic energies do not match the selected answer (Option C).\n"
        reason += f"Calculated KE of muon: {ke_mu:.2f} MeV (Expected: {expected_ke_muon} MeV)\n"
        reason += f"Calculated KE of neutrino: {ke_nu:.2f} MeV (Expected: {expected_ke_neutrino} MeV)"
        return reason

# Run the check
result = check_physics_decay_answer()
print(result)