import math

def check_correctness():
    """
    This function checks the correctness of the LLM's answer by:
    1. Calculating the theoretical kinetic energies from first principles.
    2. Verifying that the chosen option's values match the theoretical values.
    3. Explicitly checking that the chosen option satisfies both conservation of energy and momentum.
    """
    # --- Given Constants ---
    # Rest mass energies in MeV
    m_pi_c2 = 139.6
    m_mu_c2 = 105.7

    # --- LLM's Answer ---
    llm_answer_option = "D"
    options = {
        "A": (2.84, 26.8),
        "B": (3.52, 20.8),
        "C": (7.2, 32.8),
        "D": (4.12, 29.8)
    }

    # --- Theoretical Calculation ---
    # From conservation of energy and momentum, the kinetic energy of the neutrino is:
    # KE_nu = pc = ((m_pi*c^2)^2 - (m_mu*c^2)^2) / (2 * m_pi*c^2)
    ke_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)

    # The total kinetic energy released (Q-value) is the difference in rest mass energies.
    q_value = m_pi_c2 - m_mu_c2
    
    # The kinetic energy of the muon is the remaining energy.
    ke_mu_calc = q_value - ke_nu_calc

    # --- Verification ---
    if llm_answer_option not in options:
        return f"Incorrect. The chosen option '{llm_answer_option}' is not a valid option."

    # Get the values from the chosen option.
    # The lighter particle (neutrino) gets more kinetic energy.
    val1, val2 = options[llm_answer_option]
    ke_mu_llm = min(val1, val2)
    ke_nu_llm = max(val1, val2)
    
    # Set a tolerance for floating point comparisons (e.g., 0.05 MeV)
    tolerance = 0.05

    # Constraint 1: Check if the values in the option match the direct calculation.
    if not (math.isclose(ke_mu_llm, ke_mu_calc, abs_tol=tolerance) and 
            math.isclose(ke_nu_llm, ke_nu_calc, abs_tol=tolerance)):
        return (f"Incorrect. The kinetic energies in option {llm_answer_option} ({ke_mu_llm} MeV, {ke_nu_llm} MeV) "
                f"do not match the calculated theoretical values ({ke_mu_calc:.2f} MeV, {ke_nu_calc:.2f} MeV).")

    # Constraint 2: Check for Conservation of Energy (Q-value).
    # The sum of the kinetic energies must equal the total energy released.
    total_ke_llm = ke_mu_llm + ke_nu_llm
    if not math.isclose(total_ke_llm, q_value, abs_tol=tolerance):
        return (f"Incorrect. The sum of kinetic energies in option {llm_answer_option} ({total_ke_llm:.2f} MeV) "
                f"does not satisfy energy conservation. The total available kinetic energy (Q-value) is {q_value:.2f} MeV.")

    # Constraint 3: Check for Conservation of Momentum.
    # The magnitudes of the momenta of the two particles must be equal.
    # For the massless neutrino: pc_nu = KE_nu
    pc_nu = ke_nu_llm
    # For the massive muon: (pc_mu)^2 = E_mu^2 - (m_mu*c^2)^2 = (KE_mu + m_mu*c^2)^2 - (m_mu*c^2)^2
    e_mu_llm = ke_mu_llm + m_mu_c2
    pc_mu_squared = e_mu_llm**2 - m_mu_c2**2
    if pc_mu_squared < 0:
        return f"Incorrect. The values in option {llm_answer_option} lead to an invalid momentum calculation for the muon."
    pc_mu = math.sqrt(pc_mu_squared)

    if not math.isclose(pc_nu, pc_mu, abs_tol=tolerance):
        return (f"Incorrect. The momenta of the particles in option {llm_answer_option} are not equal, violating momentum conservation. "
                f"Momentum of neutrino (pc) is {pc_nu:.2f} MeV. "
                f"Calculated momentum of muon (pc) is {pc_mu:.2f} MeV.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_correctness()
print(result)