import math

def check_particle_decay_ke():
    """
    Checks the correctness of the kinetic energies of product particles in the decay
    Pi(+) -> mu(+) + nu, where Pi(+) is stationary.

    The function calculates the expected kinetic energies based on conservation laws
    and compares them to the values given in the proposed answer (Option B).
    """
    # --- Given Constants from the Question ---
    # Rest mass energies in MeV
    m_pi_c2 = 139.6
    m_mu_c2 = 105.7
    # The muon neutrino is assumed to be massless (m_nu_c2 = 0)
    # which is a standard approximation in such problems.

    # --- Answer to be Checked (Option B) ---
    # The problem asks for KE of mu(+) and nu. The heavier particle (muon)
    # will have lower kinetic energy.
    expected_ke_mu = 4.12  # MeV
    expected_ke_nu = 29.8  # MeV

    # --- Calculation based on Physics Principles ---
    # From conservation of energy and momentum, we can derive the momentum of the products.
    # E_pi = E_mu + E_nu
    # m_pi*c^2 = sqrt((p*c)^2 + (m_mu*c^2)^2) + p*c
    # Solving for p*c gives:
    # p*c = ( (m_pi*c^2)^2 - (m_mu*c^2)^2 ) / ( 2 * m_pi*c^2 )
    
    try:
        # Calculate the momentum term (p*c)
        pc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

        # For a massless neutrino, its kinetic energy is equal to its total energy, which is p*c.
        calculated_ke_nu = pc

        # For the massive muon, its total energy is E_mu = sqrt((p*c)^2 + (m_mu*c^2)^2)
        # Its kinetic energy is KE_mu = E_mu - m_mu*c^2
        e_mu = math.sqrt(pc**2 + m_mu_c2**2)
        calculated_ke_mu = e_mu - m_mu_c2
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification ---
    # We use a tolerance to account for potential rounding in the question's options.
    # A tolerance of 0.02 MeV is reasonable given the precision of the options.
    tolerance = 0.02

    # Check 1: Kinetic Energy of the muon (mu+)
    if not math.isclose(calculated_ke_mu, expected_ke_mu, abs_tol=tolerance):
        return (f"Incorrect kinetic energy for the muon (mu+).\n"
                f"Constraint: Conservation of energy and momentum.\n"
                f"Reason: The calculated KE for the muon is {calculated_ke_mu:.3f} MeV, "
                f"which does not match the expected value of {expected_ke_mu} MeV within a tolerance of {tolerance} MeV.")

    # Check 2: Kinetic Energy of the neutrino (nu)
    if not math.isclose(calculated_ke_nu, expected_ke_nu, abs_tol=tolerance):
        return (f"Incorrect kinetic energy for the neutrino (nu).\n"
                f"Constraint: Conservation of energy and momentum.\n"
                f"Reason: The calculated KE for the neutrino is {calculated_ke_nu:.3f} MeV, "
                f"which does not match the expected value of {expected_ke_nu} MeV within a tolerance of {tolerance} MeV. "
                f"The provided answer likely rounded the more precise value of ~29.78 MeV.")

    # Check 3: Sanity check on total energy conservation
    # The total kinetic energy released must equal the difference in rest mass energy.
    q_value = m_pi_c2 - m_mu_c2
    total_ke_from_answer = expected_ke_mu + expected_ke_nu
    if not math.isclose(q_value, total_ke_from_answer, abs_tol=tolerance):
         return (f"The answer violates energy conservation.\n"
                 f"Constraint: Total KE must equal the Q-value of the reaction.\n"
                 f"Reason: The Q-value (m_pi - m_mu)c^2 is {q_value:.2f} MeV. "
                 f"The sum of kinetic energies in the answer is {total_ke_from_answer:.2f} MeV. These should be equal.")


    return "Correct"

# Execute the check and print the result
result = check_particle_decay_ke()
print(result)