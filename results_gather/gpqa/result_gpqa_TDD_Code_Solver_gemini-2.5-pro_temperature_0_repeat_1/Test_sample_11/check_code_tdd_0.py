import math

def check_pion_decay_answer():
    """
    This function verifies the kinetic energies of the products in the decay
    Pi(+) -> mu(+) + nu, where the pion is stationary.

    The function recalculates the kinetic energies based on the given rest masses
    and compares them to the values provided in the proposed answer (Option B).
    """

    # --- Given Information & Constraints ---
    # Rest mass of Pi(+) in MeV/c^2
    m_pi = 139.6
    # Rest mass of mu(+) in MeV/c^2
    m_mu = 105.7
    # The neutrino (nu) is assumed to be massless, a standard assumption in such problems.
    m_nu = 0.0

    # --- The Answer to be Checked (Option B) ---
    # The question provides options with two values. By convention and calculation,
    # the smaller value is the muon's KE and the larger is the neutrino's KE.
    # KE_mu(+) = 4.12 MeV, KE_nu = 29.8 MeV
    expected_ke_mu = 4.12
    expected_ke_nu = 29.8

    # --- Physics Calculation ---
    # We use natural units where the speed of light, c = 1.
    # This means mass, momentum, and energy are all expressed in units of MeV.

    # 1. Conservation of Energy:
    # E_initial = E_final
    # Since the pion is at rest, its total energy is its rest mass energy.
    # m_pi = E_mu + E_nu

    # 2. Conservation of Momentum:
    # p_initial = p_final
    # Since the pion is at rest, its initial momentum is 0.
    # 0 = p_mu + p_nu  =>  p_mu = -p_nu
    # The magnitudes of their momenta are equal: |p_mu| = |p_nu| = p

    # 3. Relativistic Energy-Momentum Relation: E^2 = (p*c)^2 + (m*c^2)^2
    # For the muon: E_mu = sqrt(p^2 + m_mu^2)
    # For the massless neutrino: E_nu = p

    # 4. Solving for momentum 'p':
    # Substitute E_mu and E_nu into the energy conservation equation:
    # m_pi = sqrt(p^2 + m_mu^2) + p
    # Rearranging and solving for p gives the well-known formula for two-body decay:
    # p = (m_pi^2 - m_mu^2) / (2 * m_pi)

    try:
        if m_pi <= m_mu + m_nu:
            return "The decay is not energetically possible as the parent mass is not greater than the product masses."

        # Calculate the momentum 'p' of the decay products.
        p = (m_pi**2 - m_mu**2) / (2 * m_pi)

        # The kinetic energy of the massless neutrino is its total energy, which is 'p'.
        calculated_ke_nu = p

        # The total energy of the muon is calculated from the energy-momentum relation.
        E_mu = math.sqrt(p**2 + m_mu**2)
        # The kinetic energy of the muon is its total energy minus its rest mass energy.
        calculated_ke_mu = E_mu - m_mu

    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification Step ---
    # Compare the calculated values with the expected answer from Option B.
    # We use math.isclose() to handle potential floating-point inaccuracies.
    # The absolute tolerance is set to 0.01 to account for the rounding in the answer.
    is_mu_correct = math.isclose(calculated_ke_mu, expected_ke_mu, abs_tol=0.01)
    is_nu_correct = math.isclose(calculated_ke_nu, expected_ke_nu, abs_tol=0.01)

    if is_mu_correct and is_nu_correct:
        return "Correct"
    else:
        error_msg = "The answer is incorrect.\n"
        error_msg += f"Calculated values (KE_mu, KE_nu): ({calculated_ke_mu:.3f} MeV, {calculated_ke_nu:.3f} MeV)\n"
        error_msg += f"Provided answer (KE_mu, KE_nu): ({expected_ke_mu:.2f} MeV, {expected_ke_nu:.2f} MeV)\n"
        if not is_mu_correct:
            error_msg += f"Constraint failed: The muon's kinetic energy does not match.\n"
        if not is_nu_correct:
            error_msg += f"Constraint failed: The neutrino's kinetic energy does not match."
        return error_msg

# Execute the check and print the result.
result = check_pion_decay_answer()
print(result)