import math

def check_particle_decay_kinetics():
    """
    Checks the correctness of the kinetic energies of product particles in the decay
    Pi(+) -> mu(+) + nu, where Pi(+) is stationary.
    """
    # --- Given Constants from the Question ---
    # Rest mass energy of the positive pion (Pi+) in MeV
    m_pi_c2 = 139.6
    # Rest mass energy of the positive muon (mu+) in MeV
    m_mu_c2 = 105.7
    # Rest mass energy of the neutrino (nu) is assumed to be 0 MeV, a standard approximation.
    m_nu_c2 = 0

    # --- Answer to be Checked (Option B) ---
    # Kinetic energy of the muon in MeV
    answer_KE_mu = 4.12
    # Kinetic energy of the neutrino in MeV
    answer_KE_nu = 29.8

    # --- Physics Calculations ---
    # In the rest frame of the pion, the initial momentum is 0. By conservation of momentum,
    # the muon and neutrino must have equal and opposite momenta. Let the magnitude be 'p'.
    #
    # By conservation of energy:
    # E_pi = E_mu + E_nu
    # m_pi_c2 = sqrt((p*c)^2 + (m_mu_c2)^2) + sqrt((p*c)^2 + (m_nu_c2)^2)
    # Since m_nu_c2 = 0, this simplifies to:
    # m_pi_c2 = sqrt((p*c)^2 + (m_mu_c2)^2) + p*c
    #
    # We can solve this equation for the momentum term, p*c:
    # p*c = ( (m_pi_c2)^2 - (m_mu_c2)^2 ) / ( 2 * m_pi_c2 )

    try:
        # Calculate the momentum-energy term (pc) for the products
        pc_numerator = m_pi_c2**2 - m_mu_c2**2
        pc_denominator = 2 * m_pi_c2
        pc = pc_numerator / pc_denominator

        # --- Calculate Kinetic Energies from Momentum ---

        # For the massless neutrino, Kinetic Energy = Total Energy = pc
        calculated_KE_nu = pc

        # For the massive muon, Kinetic Energy = Total Energy - Rest Mass Energy
        # Total Energy E_mu = sqrt((pc)^2 + (m_mu_c2)^2)
        total_energy_mu = math.sqrt(pc**2 + m_mu_c2**2)
        calculated_KE_mu = total_energy_mu - m_mu_c2

    except Exception as e:
        return f"An error occurred during the physics calculation: {e}"

    # --- Verification Step ---
    # Compare the calculated values with the provided answer, allowing for a small
    # tolerance due to rounding in the problem's options.
    # An absolute tolerance of 0.02 is appropriate given the answer's precision.
    tolerance = 0.02

    is_muon_ke_correct = math.isclose(calculated_KE_mu, answer_KE_mu, abs_tol=tolerance)
    is_neutrino_ke_correct = math.isclose(calculated_KE_nu, answer_KE_nu, abs_tol=tolerance)

    if is_muon_ke_correct and is_neutrino_ke_correct:
        return "Correct"
    else:
        error_message = "The provided answer is incorrect.\n"
        if not is_muon_ke_correct:
            error_message += (f"Constraint check failed for muon kinetic energy (KE_mu).\n"
                              f"  - Calculated KE_mu: {calculated_KE_mu:.3f} MeV\n"
                              f"  - Provided answer KE_mu: {answer_KE_mu} MeV\n")
        if not is_neutrino_ke_correct:
            error_message += (f"Constraint check failed for neutrino kinetic energy (KE_nu).\n"
                              f"  - Calculated KE_nu: {calculated_KE_nu:.3f} MeV\n"
                              f"  - Provided answer KE_nu: {answer_KE_nu} MeV\n")
        
        # Final sanity check: The sum of kinetic energies must equal the Q-value of the decay.
        q_value = m_pi_c2 - m_mu_c2 - m_nu_c2
        total_calculated_ke = calculated_KE_mu + calculated_KE_nu
        if not math.isclose(q_value, total_calculated_ke, rel_tol=1e-9):
             error_message += (f"Energy conservation check failed: The sum of calculated kinetic energies "
                               f"({total_calculated_ke:.3f} MeV) does not equal the reaction's Q-value ({q_value:.1f} MeV).")

        return error_message.strip()

# Execute the check
result = check_particle_decay_kinetics()
print(result)