import math

def check_particle_decay_ke():
    """
    Checks the correctness of the answer for the kinetic energies in pion decay.

    The decay is Pi(+) -> mu(+) + nu, with the Pi(+) being stationary.
    Rest mass of Pi(+) (m_pi_c2) = 139.6 MeV
    Rest mass of mu(+) (m_mu_c2) = 105.7 MeV
    The neutrino (nu) is assumed to be massless.

    The final answer to be checked is C, which corresponds to:
    KE_mu = 4.12 MeV, KE_nu = 29.8 MeV
    """

    # 1. Define the physical constants from the question
    m_pi_c2 = 139.6  # Rest energy of the pion in MeV
    m_mu_c2 = 105.7  # Rest energy of the muon in MeV

    # 2. Calculate the theoretical kinetic energies from first principles
    # From conservation of energy and momentum, the kinetic energy of the neutrino (KE_nu)
    # can be calculated as: KE_nu = pc = ((m_pi_c2^2) - (m_mu_c2^2)) / (2 * m_pi_c2)
    try:
        ke_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation error: Division by zero. The rest mass of the pion cannot be zero."

    # The total kinetic energy released (Q-value) is the difference in rest masses.
    # Q = KE_mu + KE_nu = m_pi_c2 - m_mu_c2
    q_value = m_pi_c2 - m_mu_c2
    
    # The kinetic energy of the muon is the remaining kinetic energy.
    ke_mu_calc = q_value - ke_nu_calc

    # 3. Define the values from the proposed answer (Option C)
    # The heavier particle (muon) will have lower kinetic energy.
    answer_ke_mu = 4.12
    answer_ke_nu = 29.8

    # 4. Check if the calculated values match the answer's values
    # We use a tolerance to account for potential rounding in the options.
    tolerance = 0.01 

    # Check 1: Compare calculated KE_mu with the answer's KE_mu
    if not math.isclose(ke_mu_calc, answer_ke_mu, rel_tol=tolerance):
        return (f"Incorrect. The calculated kinetic energy of the muon is {ke_mu_calc:.2f} MeV, "
                f"but the answer C provides {answer_ke_mu} MeV. The values do not match.")

    # Check 2: Compare calculated KE_nu with the answer's KE_nu
    if not math.isclose(ke_nu_calc, answer_ke_nu, rel_tol=tolerance):
        return (f"Incorrect. The calculated kinetic energy of the neutrino is {ke_nu_calc:.2f} MeV, "
                f"but the answer C provides {answer_ke_nu} MeV. The values do not match.")

    # 5. Verify that the answer itself satisfies the conservation laws (as a sanity check)
    # Check conservation of energy (Q-value)
    if not math.isclose(answer_ke_mu + answer_ke_nu, q_value, rel_tol=tolerance):
        return (f"Incorrect. The sum of kinetic energies in answer C ({answer_ke_mu + answer_ke_nu:.2f} MeV) "
                f"does not match the total kinetic energy released ({q_value:.2f} MeV).")

    # Check conservation of momentum: |p_mu| must equal |p_nu|
    # For the muon: (p_mu*c)^2 = E_mu^2 - (m_mu*c^2)^2 = (KE_mu + m_mu*c^2)^2 - (m_mu*c^2)^2
    # For the massless neutrino: (p_nu*c)^2 = E_nu^2 = KE_nu^2
    p_mu_c_sq = (answer_ke_mu + m_mu_c2)**2 - m_mu_c2**2
    p_nu_c_sq = answer_ke_nu**2
    
    if not math.isclose(p_mu_c_sq, p_nu_c_sq, rel_tol=tolerance):
        return (f"Incorrect. The momenta of the particles in answer C are not conserved. "
                f"(p_mu*c)^2 is {p_mu_c_sq:.2f} MeV^2 while (p_nu*c)^2 is {p_nu_c_sq:.2f} MeV^2.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_particle_decay_ke()
print(result)