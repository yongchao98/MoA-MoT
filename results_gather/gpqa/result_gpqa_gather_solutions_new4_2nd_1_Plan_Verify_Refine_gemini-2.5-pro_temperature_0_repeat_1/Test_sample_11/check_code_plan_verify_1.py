import math

def check_answer():
    """
    Checks the correctness of the LLM's answer for the particle decay problem.
    """
    # Given values from the question
    m_pi_c2 = 139.6  # Rest energy of Pi(+) in MeV
    m_mu_c2 = 105.7  # Rest energy of mu(+) in MeV

    # Options provided in the question
    options = {
        'A': [7.2, 32.8],
        'B': [3.52, 20.8],
        'C': [4.12, 29.8],
        'D': [2.84, 26.8]
    }

    # The final answer provided by the LLM
    llm_answer_letter = 'D'

    # --- Physics Calculation ---
    # Calculate the kinetic energy of the neutrino (KE_nu)
    # KE_nu = pc = ((m_pi*c^2)^2 - (m_mu*c^2)^2) / (2 * m_pi*c^2)
    ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

    # Calculate the total kinetic energy released (Q-value)
    # KE_total = m_pi*c^2 - m_mu*c^2
    ke_total = m_pi_c2 - m_mu_c2

    # Calculate the kinetic energy of the muon (KE_mu)
    # KE_mu = KE_total - KE_nu
    ke_mu = ke_total - ke_nu

    calculated_values = sorted([ke_mu, ke_nu])

    # --- Verification ---
    # Get the values from the LLM's chosen option
    llm_answer_values = sorted(options.get(llm_answer_letter))

    if not llm_answer_values:
        return f"Invalid answer format. The selected option '{llm_answer_letter}' does not exist."

    # Compare the calculated values with the values from the chosen option
    # We use math.isclose to handle potential floating-point rounding differences.
    # A tolerance of 0.01 MeV is reasonable for this problem.
    is_correct = (math.isclose(calculated_values[0], llm_answer_values[0], abs_tol=0.01) and
                  math.isclose(calculated_values[1], llm_answer_values[1], abs_tol=0.01))

    if is_correct:
        return "Correct"
    else:
        return (f"Incorrect. The final answer is D, which corresponds to the values {options['D']} MeV. "
                f"However, the calculation based on conservation of energy and momentum yields kinetic energies of "
                f"approximately {ke_mu:.2f} MeV for the muon and {ke_nu:.2f} MeV for the neutrino. "
                f"These calculated values do not match the values in option D.")

# Run the check
result = check_answer()
print(result)