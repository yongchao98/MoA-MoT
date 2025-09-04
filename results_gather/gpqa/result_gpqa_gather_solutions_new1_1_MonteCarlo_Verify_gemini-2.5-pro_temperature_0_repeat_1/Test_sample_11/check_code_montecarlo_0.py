import math

def check_particle_decay_answer():
    """
    Checks the correctness of the kinetic energies for the decay Pi(+) -> mu(+) + nu.
    The final answer provided is D, with KE_mu = 4.12 MeV and KE_nu = 29.8 MeV.
    """
    # Given physical constants from the question
    m_pi_c2 = 139.6  # Rest energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest energy of mu(+) in MeV

    # The proposed answer values from option D
    KE_mu_proposed = 4.12  # Proposed kinetic energy of the muon in MeV
    KE_nu_proposed = 29.8  # Proposed kinetic energy of the neutrino in MeV

    # Tolerance for comparing floating-point numbers due to rounding in the options
    tolerance = 1e-2

    # --- 1. Check Conservation of Energy ---
    # The total kinetic energy released (Q-value) must equal the difference in rest energies.
    q_value_theory = m_pi_c2 - m_mu_c2
    q_value_proposed = KE_mu_proposed + KE_nu_proposed

    if not math.isclose(q_value_theory, q_value_proposed, rel_tol=tolerance):
        return (f"Incorrect. Conservation of energy is not satisfied. "
                f"The sum of the proposed kinetic energies is {q_value_proposed:.2f} MeV, "
                f"but the theoretical Q-value (m_pi*c^2 - m_mu*c^2) is {q_value_theory:.2f} MeV.")

    # --- 2. Check Conservation of Momentum ---
    # The magnitudes of the momenta of the two product particles must be equal.
    # We will compare the square of their momenta, (p*c)^2, to avoid using square roots.

    # For the massive muon, from E^2 = (pc)^2 + (mc^2)^2:
    # (p_mu*c)^2 = E_mu^2 - (m_mu*c^2)^2
    # where E_mu = KE_mu + m_mu*c^2
    E_mu_proposed = KE_mu_proposed + m_mu_c2
    p_mu_c_sq = E_mu_proposed**2 - m_mu_c2**2

    # For the massless neutrino, E_nu = p_nu*c.
    # Since its total energy is its kinetic energy, E_nu = KE_nu.
    # So, (p_nu*c)^2 = KE_nu^2
    p_nu_c_sq = KE_nu_proposed**2

    if not math.isclose(p_mu_c_sq, p_nu_c_sq, rel_tol=tolerance):
        return (f"Incorrect. Conservation of momentum is not satisfied. "
                f"The momentum magnitudes are not equal. "
                f"Calculated (p_mu*c)^2 = {p_mu_c_sq:.2f} MeV^2. "
                f"Calculated (p_nu*c)^2 = {p_nu_c_sq:.2f} MeV^2.")

    # If both fundamental constraints are satisfied, the answer is correct.
    return "Correct"

# Run the check and print the result
result = check_particle_decay_answer()
print(result)