import math

def check_answer():
    """
    Checks the correctness of the calculated kinetic energies for the pion decay.
    Pi(+) = mu(+) + nu
    """
    # Given values from the question
    m_pi_c2 = 139.6  # Rest mass energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest mass energy of mu(+) in MeV
    m_nu_c2 = 0      # Assumed rest mass energy of neutrino in MeV

    # The answer to be checked, taken from the provided LLM response.
    # The response identifies option B as correct.
    # Option B: KE_mu = 4.12 MeV, KE_nu = 29.8 MeV
    # The order in the option is not specified, so we must check both possibilities.
    # However, the muon is much heavier, so it will have the smaller kinetic energy.
    KE_mu_ans = 4.12
    KE_nu_ans = 29.8

    # --- Verification Step 1: Calculate the correct values from first principles ---
    # This is the most robust way to check. We solve the problem from scratch.
    # From conservation of energy and momentum, we can derive the energy of the neutrino.
    # E_nu = pc = ((m_pi^2 - m_mu^2) * c^4) / (2 * m_pi * c^2)
    
    try:
        # Calculate kinetic energy of the neutrino
        KE_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

        # The total kinetic energy released (Q-value) is the difference in rest mass energies
        Q_value = m_pi_c2 - m_mu_c2

        # The remaining kinetic energy must belong to the muon
        KE_mu_calc = Q_value - KE_nu_calc
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # --- Verification Step 2: Compare calculated values with the provided answer ---
    # We use a small tolerance to account for potential rounding in the answer.
    tolerance = 0.01

    # Check if the calculated muon KE matches the answer's muon KE
    if not math.isclose(KE_mu_calc, KE_mu_ans, rel_tol=tolerance):
        return (f"Incorrect. The kinetic energy of the muon (mu+) is incorrect.\n"
                f"Expected value based on conservation laws: ~{KE_mu_calc:.2f} MeV.\n"
                f"Provided answer: {KE_mu_ans} MeV.")

    # Check if the calculated neutrino KE matches the answer's neutrino KE
    if not math.isclose(KE_nu_calc, KE_nu_ans, rel_tol=tolerance):
        return (f"Incorrect. The kinetic energy of the neutrino (nu) is incorrect.\n"
                f"Expected value based on conservation laws: ~{KE_nu_calc:.2f} MeV.\n"
                f"Provided answer: {KE_nu_ans} MeV.")

    # --- Verification Step 3: Cross-check using conservation laws on the answer's values ---
    # This confirms the answer is self-consistent.

    # Constraint 1: Conservation of Energy (Total KE must equal Q-value)
    Q_value_ans = KE_mu_ans + KE_nu_ans
    if not math.isclose(Q_value_ans, Q_value, rel_tol=tolerance):
        return (f"Incorrect. The answer violates conservation of energy.\n"
                f"The sum of kinetic energies in the answer ({KE_mu_ans} + {KE_nu_ans} = {Q_value_ans:.2f} MeV) "
                f"does not equal the available energy from mass difference ({Q_value:.2f} MeV).")

    # Constraint 2: Conservation of Momentum (Momenta must be equal and opposite)
    # (p_mu * c)^2 = E_mu^2 - (m_mu * c^2)^2 = (KE_mu + m_mu*c^2)^2 - (m_mu*c^2)^2
    # (p_nu * c)^2 = E_nu^2 - (m_nu * c^2)^2 = KE_nu^2
    p_mu_c_sq = (KE_mu_ans + m_mu_c2)**2 - m_mu_c2**2
    p_nu_c_sq = KE_nu_ans**2
    
    if not math.isclose(p_mu_c_sq, p_nu_c_sq, rel_tol=tolerance):
        return (f"Incorrect. The answer violates conservation of momentum.\n"
                f"The momentum of the muon and neutrino are not equal.\n"
                f"Calculated (p_mu*c)^2 = {p_mu_c_sq:.2f} MeV^2.\n"
                f"Calculated (p_nu*c)^2 = {p_nu_c_sq:.2f} MeV^2.")

    return "Correct"

# Run the check
result = check_answer()
print(result)