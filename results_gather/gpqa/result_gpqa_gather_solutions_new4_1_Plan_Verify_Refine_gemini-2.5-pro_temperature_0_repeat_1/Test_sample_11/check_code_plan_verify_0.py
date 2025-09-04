import math

def check_particle_decay_ke():
    """
    Checks the correctness of the kinetic energies of product particles in a pion decay.

    The decay process is: Pi(+) -> mu(+) + nu, where Pi(+) is stationary.
    This function calculates the kinetic energies based on conservation of energy and momentum
    and compares them to the values given in the proposed answer.
    """
    # --- Given values from the question ---
    m_pi_c2 = 139.6  # Rest mass energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest mass energy of mu(+) in MeV
    # The neutrino (nu) is assumed to be massless.

    # --- Values from the proposed answer (Option D) ---
    # The answer states the KE of the products are 4.12 MeV and 29.8 MeV.
    # We need to check which value corresponds to which particle.
    expected_values = {4.12, 29.8}

    # --- Calculation based on Relativistic Kinematics ---
    # For a two-body decay of a stationary particle, we can derive the energies
    # of the products from conservation of energy and momentum.

    # The energy of the neutrino (E_nu) can be calculated with the formula:
    # E_nu = [ (m_pi_c2)^2 - (m_mu_c2)^2 ] / [ 2 * m_pi_c2 ]
    # Since the neutrino is massless, its total energy is its kinetic energy (KE_nu).
    try:
        calculated_ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation Error: The rest mass of the initial particle (pion) cannot be zero."

    # The total kinetic energy released (Q-value) is the difference in rest masses.
    # Q = m_pi_c2 - m_mu_c2
    # This energy is shared between the two product particles: Q = KE_mu + KE_nu
    # Therefore, the kinetic energy of the muon (KE_mu) is:
    calculated_ke_mu = (m_pi_c2 - m_mu_c2) - calculated_ke_nu

    # --- Verification ---
    # We round the calculated values to two decimal places to match the precision of the answer.
    calculated_values = {round(calculated_ke_mu, 2), round(calculated_ke_nu, 2)}

    # Check if the set of calculated values matches the set of expected values.
    # Using sets makes the comparison independent of the order of the values.
    if calculated_values == expected_values:
        return "Correct"
    else:
        # If the values don't match, provide a detailed reason.
        reason = (
            f"Incorrect. The calculated kinetic energies do not match the answer's values.\n"
            f"The answer provides the pair of kinetic energies: {expected_values} MeV.\n"
            f"The calculation based on conservation laws yields:\n"
            f"  - KE of muon (mu+): {calculated_ke_mu:.2f} MeV\n"
            f"  - KE of neutrino (nu): {calculated_ke_nu:.2f} MeV\n"
            f"The calculated pair of values is {calculated_values} MeV, which does not match the answer."
        )
        return reason

# Execute the check and print the result.
print(check_particle_decay_ke())