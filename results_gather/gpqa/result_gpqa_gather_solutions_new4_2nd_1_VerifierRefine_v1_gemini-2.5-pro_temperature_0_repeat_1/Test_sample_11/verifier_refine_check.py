import math

def check_correctness():
    """
    This function checks the correctness of the provided answer for the particle physics decay problem.
    It calculates the kinetic energies of the product particles based on relativistic conservation laws
    and compares them to the values given in the selected answer option.
    """

    # --- Problem Constraints and Given Values ---
    # Rest mass of Pi(+) in MeV/c^2 (equivalent to rest energy in MeV)
    m_pi_c2 = 139.6
    # Rest mass of mu(+) in MeV/c^2 (equivalent to rest energy in MeV)
    m_mu_c2 = 105.7
    # The neutrino is assumed to be massless.
    m_nu_c2 = 0

    # --- The Answer to be Checked ---
    # The final answer provided is 'D', which corresponds to the values:
    # KE_mu = 4.12 MeV
    # KE_nu = 29.8 MeV
    # Note: In a two-body decay, the lighter particle (neutrino) gets more kinetic energy.
    expected_ke_mu = 4.12
    expected_ke_nu = 29.8

    # --- Physics Calculation ---
    # The calculation is based on conservation of energy and momentum in a relativistic framework.
    # For a two-body decay of a stationary particle (M) into two particles (m1, m2),
    # the energy of particle 1 (E1) can be calculated as:
    # E1 = (M^2 + m1^2 - m2^2) / (2 * M)
    # Here, M = m_pi_c2, m1 = m_mu_c2, m2 = m_nu_c2 = 0
    # The total energy of the muon (E_mu) is:
    e_mu_calc = (m_pi_c2**2 + m_mu_c2**2 - m_nu_c2**2) / (2 * m_pi_c2)

    # The kinetic energy of the muon is its total energy minus its rest energy.
    ke_mu_calc = e_mu_calc - m_mu_c2

    # The total kinetic energy released (Q-value) is the difference in rest masses.
    ke_total = m_pi_c2 - m_mu_c2

    # This total kinetic energy is shared between the two particles.
    # KE_total = KE_mu + KE_nu
    # Therefore, the kinetic energy of the neutrino is:
    ke_nu_calc = ke_total - ke_mu_calc

    # --- Verification Step ---
    # Check if the calculated values match the expected values from the answer.
    # A tolerance is used for floating-point comparisons.
    tolerance = 1e-2  # Corresponds to 0.01 MeV

    muon_ke_match = math.isclose(ke_mu_calc, expected_ke_mu, rel_tol=tolerance)
    neutrino_ke_match = math.isclose(ke_nu_calc, expected_ke_nu, rel_tol=tolerance)

    if muon_ke_match and neutrino_ke_match:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated kinetic energies do not match the values in the selected answer (D).\n"
        error_message += f"Constraint Check:\n"
        error_message += f"1. Calculated KE of muon (mu+): {ke_mu_calc:.2f} MeV. Expected: {expected_ke_mu} MeV. Match: {muon_ke_match}\n"
        error_message += f"2. Calculated KE of neutrino (nu): {ke_nu_calc:.2f} MeV. Expected: {expected_ke_nu} MeV. Match: {neutrino_ke_match}\n"
        error_message += f"The physics calculation yields KE_mu ≈ 4.12 MeV and KE_nu ≈ 29.8 MeV, which should match the answer."
        return error_message

# The final answer from the LLM is <<<D>>>, which corresponds to the values 4.12 MeV and 29.8 MeV.
# The code will now check if these values are correct based on the physics principles.
result = check_correctness()
print(result)