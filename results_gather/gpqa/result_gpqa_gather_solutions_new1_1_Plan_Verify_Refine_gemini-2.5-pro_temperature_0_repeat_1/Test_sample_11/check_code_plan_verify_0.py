import math

def check_particle_decay_ke():
    """
    Calculates the kinetic energies of the product particles in pion decay
    and checks them against the provided answer.
    """
    # Given values from the question
    m_pi_c2 = 139.6  # Rest mass energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest mass energy of mu(+) in MeV

    # The final answer from the LLM analysis is <<<B>>>.
    # Option B is: 4.12 MeV, 29.8 MeV
    # By convention and calculation, the smaller energy belongs to the massive muon.
    expected_ke_mu = 4.12
    expected_ke_nu = 29.8

    # --- Calculation from First Principles ---

    # 1. Calculate the kinetic energy of the massless neutrino (KE_nu).
    # This is derived from conservation of energy and momentum.
    # KE_nu = ( (m_pi_c2)^2 - (m_mu_c2)^2 ) / ( 2 * m_pi_c2 )
    try:
        numerator = m_pi_c2**2 - m_mu_c2**2
        denominator = 2 * m_pi_c2
        calculated_ke_nu = numerator / denominator
    except ZeroDivisionError:
        return "Error: Division by zero. The rest mass of the pion cannot be zero."

    # 2. Calculate the total kinetic energy released (Q-value).
    q_value = m_pi_c2 - m_mu_c2

    # 3. Calculate the kinetic energy of the muon (KE_mu).
    # The Q-value is shared between the two particles: Q = KE_mu + KE_nu
    calculated_ke_mu = q_value - calculated_ke_nu

    # --- Verification ---

    # Set a tolerance for floating-point comparison, e.g., 1% relative tolerance.
    tolerance = 0.01

    # Check if the calculated values match the expected answer within the tolerance.
    muon_ke_matches = math.isclose(calculated_ke_mu, expected_ke_mu, rel_tol=tolerance)
    neutrino_ke_matches = math.isclose(calculated_ke_nu, expected_ke_nu, rel_tol=tolerance)

    if muon_ke_matches and neutrino_ke_matches:
        return "Correct"
    else:
        reason = "Incorrect. The provided answer does not match the values calculated from physical principles.\n"
        reason += f"Calculated KE for muon: {calculated_ke_mu:.3f} MeV. Expected: {expected_ke_mu} MeV.\n"
        reason += f"Calculated KE for neutrino: {calculated_ke_nu:.3f} MeV. Expected: {expected_ke_nu} MeV.\n"
        if not muon_ke_matches:
            reason += "The kinetic energy for the muon is incorrect.\n"
        if not neutrino_ke_matches:
            reason += "The kinetic energy for the neutrino is incorrect.\n"
        return reason

# Run the check
result = check_particle_decay_ke()
print(result)