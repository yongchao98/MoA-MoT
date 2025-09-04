import math

def check_particle_decay_answer():
    """
    Checks the correctness of the kinetic energies of product particles in the decay
    Pi(+) -> mu(+) + nu, where Pi(+) is stationary.
    """
    # --- Given Constants ---
    # Rest mass energy of the positive pion (Pi+) in MeV
    m_pi_c2 = 139.6
    # Rest mass energy of the positive muon (mu+) in MeV
    m_mu_c2 = 105.7
    # The neutrino (nu) is considered massless, so its rest mass energy is 0.

    # --- Proposed Answer (Option B) ---
    # Kinetic energy of the muon in MeV
    KE_mu_proposed = 4.12
    # Kinetic energy of the neutrino in MeV
    KE_nu_proposed = 29.8

    # --- Set a tolerance for floating-point comparisons ---
    # The values in the question are given to one decimal place, so a small
    # tolerance is needed to account for potential rounding.
    # An absolute tolerance of 0.1 MeV is reasonable.
    tolerance = 0.1

    # --- Verification 1: Conservation of Energy ---
    # The total kinetic energy released (Q-value) is the difference in rest mass energy.
    q_value = m_pi_c2 - m_mu_c2
    # The total kinetic energy from the proposed answer.
    total_ke_proposed = KE_mu_proposed + KE_nu_proposed

    if not math.isclose(q_value, total_ke_proposed, abs_tol=tolerance):
        return (f"Incorrect: Energy conservation is not satisfied. "
                f"The total kinetic energy released should be the mass difference: "
                f"{m_pi_c2} - {m_mu_c2} = {q_value:.2f} MeV. "
                f"The proposed answer's total kinetic energy is {KE_mu_proposed} + {KE_nu_proposed} = {total_ke_proposed:.2f} MeV.")

    # --- Verification 2: Conservation of Momentum ---
    # For a stationary decay, the momentum magnitudes of the two products must be equal.
    # We use the relativistic energy-momentum relation: E^2 = (pc)^2 + (mc^2)^2
    # which can be rearranged to find the momentum: pc = sqrt(E^2 - (mc^2)^2)

    # Calculate the muon's momentum (pc_mu)
    # Total energy of muon E_mu = KE_mu + m_mu*c^2
    E_mu = KE_mu_proposed + m_mu_c2
    # (pc_mu)^2 = E_mu^2 - (m_mu*c^2)^2
    pc_mu_squared = E_mu**2 - m_mu_c2**2
    pc_mu = math.sqrt(pc_mu_squared)

    # Calculate the neutrino's momentum (pc_nu)
    # For a massless neutrino, E_nu = pc_nu. Its energy is purely kinetic.
    pc_nu = KE_nu_proposed

    if not math.isclose(pc_mu, pc_nu, abs_tol=tolerance):
        return (f"Incorrect: Momentum conservation is not satisfied. "
                f"The magnitude of the muon's momentum (p_mu*c) is calculated to be {pc_mu:.2f} MeV. "
                f"The magnitude of the neutrino's momentum (p_nu*c) is {pc_nu:.2f} MeV. "
                f"These values should be equal for momentum to be conserved.")

    # If both conservation laws are satisfied, the answer is correct.
    return "Correct"

# Execute the check and print the result.
result = check_particle_decay_answer()
print(result)