import math

def check_pion_decay_kinetics():
    """
    This function verifies the kinetic energies of the products in the decay of a stationary pion.
    The decay is Pi(+) -> mu(+) + nu.
    The calculation is based on the principles of conservation of energy and momentum.
    """
    # --- Given values from the question ---
    # Rest energies are given in MeV.
    m_pi_c2 = 139.6  # Rest energy of the stationary positive pion (Pi+)
    m_mu_c2 = 105.7  # Rest energy of the positive muon (mu+)
    m_nu_c2 = 0.0    # The rest energy of the neutrino is assumed to be zero, a standard practice.

    # --- Answer to be checked (from option A) ---
    # The question asks for KE of product particles, mu(+) and nu.
    # Option A provides the values in the order: KE_mu, KE_nu.
    expected_KE_mu = 4.12  # MeV
    expected_KE_nu = 29.8  # MeV

    # --- Physics Calculation ---
    # From conservation of energy and momentum, one can derive the formula for the
    # momentum (times c) of either product particle.
    # p*c = [ (m_pi*c^2)^2 - (m_mu*c^2)^2 ] / [ 2 * m_pi*c^2 ]
    try:
        pc_numerator = m_pi_c2**2 - m_mu_c2**2
        pc_denominator = 2 * m_pi_c2
        pc = pc_numerator / pc_denominator
    except ZeroDivisionError:
        return "Error: The rest mass of the pion cannot be zero."

    # The kinetic energy of the massless neutrino is equal to its total energy, which is p*c.
    calculated_KE_nu = pc

    # The kinetic energy of the muon can be found from the conservation of energy:
    # Initial Energy = Final Energy
    # m_pi_c2 = E_mu + E_nu
    # m_pi_c2 = (KE_mu + m_mu_c2) + KE_nu
    # KE_mu = m_pi_c2 - m_mu_c2 - KE_nu
    calculated_KE_mu = m_pi_c2 - m_mu_c2 - calculated_KE_nu

    # --- Verification ---
    # Check if the calculated values match the expected values from the answer.
    # We use math.isclose() to handle potential floating-point inaccuracies.
    # The tolerance is set based on the precision of the provided answer (2 decimal places).
    muon_ke_correct = math.isclose(calculated_KE_mu, expected_KE_mu, abs_tol=0.01)
    neutrino_ke_correct = math.isclose(calculated_KE_nu, expected_KE_nu, abs_tol=0.01)

    if muon_ke_correct and neutrino_ke_correct:
        return "Correct"
    else:
        error_message = "The answer is incorrect.\n"
        if not muon_ke_correct:
            error_message += (f"Constraint not satisfied: Kinetic energy of the muon (mu+).\n"
                              f"Calculated value: {calculated_KE_mu:.3f} MeV, but the answer provides {expected_KE_mu} MeV.\n")
        if not neutrino_ke_correct:
            error_message += (f"Constraint not satisfied: Kinetic energy of the neutrino (nu).\n"
                              f"Calculated value: {calculated_KE_nu:.3f} MeV, but the answer provides {expected_KE_nu} MeV.\n")
        return error_message.strip()

# Execute the check and print the result.
result = check_pion_decay_kinetics()
print(result)