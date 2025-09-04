import math

def check_correctness():
    """
    Checks the correctness of the answer for the pion decay problem.

    The problem involves a stationary pion decaying into a muon and a neutrino.
    Pi(+) -> mu(+) + nu

    We use conservation of energy and momentum to calculate the kinetic energies
    of the product particles.
    """

    # Given rest mass energies from the question
    m_pi_c2 = 139.6  # MeV
    m_mu_c2 = 105.7  # MeV
    # The neutrino is assumed to be massless, m_nu_c2 = 0

    # --- Calculation based on physics principles ---

    # From conservation of energy and momentum, we can derive the formula for the
    # momentum-energy (pc) of the products, which is also the kinetic energy
    # of the massless neutrino (KE_nu).
    # KE_nu = pc = ((m_pi*c^2)^2 - (m_mu*c^2)^2) / (2 * m_pi*c^2)
    try:
        calculated_ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Error: The rest mass of the pion cannot be zero."

    # The total kinetic energy released (Q-value) is the difference in rest masses.
    # Q = m_pi*c^2 - m_mu*c^2
    # This energy is shared between the products: Q = KE_mu + KE_nu
    q_value = m_pi_c2 - m_mu_c2
    calculated_ke_mu = q_value - calculated_ke_nu

    # --- Verification against the provided answer ---

    # The final answer from the LLM is 'B', which corresponds to the option:
    # B) 4.12 MeV, 29.8 MeV
    # By convention and physics (the lighter particle gets more energy),
    # the smaller value is for the muon and the larger for the neutrino.
    answer_ke_mu = 4.12
    answer_ke_nu = 29.8

    # We use math.isclose() for a safe floating-point comparison, allowing for
    # small rounding differences. A relative tolerance of 1% is generous enough.
    is_muon_ke_correct = math.isclose(calculated_ke_mu, answer_ke_mu, rel_tol=0.01)
    is_neutrino_ke_correct = math.isclose(calculated_ke_nu, answer_ke_nu, rel_tol=0.01)

    if is_muon_ke_correct and is_neutrino_ke_correct:
        return "Correct"
    else:
        error_message = "Incorrect. The calculated kinetic energies do not match the values in the selected answer.\n"
        error_message += f"Calculated KE for muon (mu+): {calculated_ke_mu:.4f} MeV. Answer's value: {answer_ke_mu} MeV.\n"
        error_message += f"Calculated KE for neutrino (nu): {calculated_ke_nu:.4f} MeV. Answer's value: {answer_ke_nu} MeV.\n"
        if not is_muon_ke_correct:
            error_message += "The kinetic energy for the muon is incorrect.\n"
        if not is_neutrino_ke_correct:
            error_message += "The kinetic energy for the neutrino is incorrect.\n"
        return error_message

# Execute the check and print the result
result = check_correctness()
print(result)