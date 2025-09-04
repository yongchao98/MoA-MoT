import math

def check_pion_decay_answer():
    """
    This function verifies the kinetic energies of the products in the decay
    Pi(+) -> mu(+) + nu, where the initial Pi(+) is stationary.

    It checks two physical constraints:
    1. Conservation of Energy: The sum of the product KEs must equal the Q-value.
    2. Conservation of Momentum: The product momenta must be equal in magnitude.
    """
    # --- Given Constants (in MeV and MeV/c^2) ---
    m_pi = 139.6  # Rest mass of Pi(+)
    m_mu = 105.7  # Rest mass of mu(+)
    # The neutrino (nu) is assumed to be massless.

    # --- The proposed answer to check ---
    # This corresponds to option A
    ke_mu_proposed = 4.12  # Proposed KE of the muon (MeV)
    ke_nu_proposed = 29.8  # Proposed KE of the neutrino (MeV)

    # --- Verification ---

    # 1. Check Conservation of Energy
    # The total kinetic energy released (Q-value) is the difference in rest mass.
    q_value = m_pi - m_mu
    total_ke_proposed = ke_mu_proposed + ke_nu_proposed

    # Use a tolerance for floating-point comparison (e.g., 1%)
    if not math.isclose(total_ke_proposed, q_value, rel_tol=0.01):
        return (f"Incorrect: The answer violates the conservation of energy.\n"
                f"The sum of the proposed kinetic energies is {total_ke_proposed:.2f} MeV.\n"
                f"The expected total kinetic energy (Q-value) is m_pi - m_mu = {m_pi} - {m_mu} = {q_value:.2f} MeV.")

    # 2. Check Conservation of Momentum
    # Since the initial momentum is zero, the final momenta must be equal in magnitude.
    # |p_mu| = |p_nu|  =>  (p_mu*c)^2 = (p_nu*c)^2
    #
    # We use the relativistic energy-momentum relation: E^2 = (pc)^2 + (mc^2)^2
    # which gives (pc)^2 = E^2 - (mc^2)^2.
    # We work in units where c=1.

    # For the massive muon (mu+):
    # E_mu = KE_mu + m_mu
    p_mu_squared = (ke_mu_proposed + m_mu)**2 - m_mu**2

    # For the massless neutrino (nu):
    # E_nu = KE_nu and m_nu = 0, so (p_nu)^2 = E_nu^2 = KE_nu^2
    p_nu_squared = ke_nu_proposed**2

    # Compare the squared momenta using a tolerance
    if not math.isclose(p_mu_squared, p_nu_squared, rel_tol=0.01):
        return (f"Incorrect: The answer violates the conservation of momentum.\n"
                f"The magnitudes of the product momenta are not equal.\n"
                f"Calculated (p_mu*c)^2 = {p_mu_squared:.2f} (MeV)^2\n"
                f"Calculated (p_nu*c)^2 = {p_nu_squared:.2f} (MeV)^2\n"
                f"These values should be equal.")

    # If both checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_pion_decay_answer()
print(result)