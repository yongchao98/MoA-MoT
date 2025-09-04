import math

def check_correctness():
    """
    Checks the correctness of the answer for the particle decay problem.

    The problem asks for the kinetic energies (KE) of the product particles
    in the decay: Pi(+) = mu(+) + nu, where Pi(+) is stationary.

    Given:
    - Rest mass of Pi(+) (m_pi_c2) = 139.6 MeV
    - Rest mass of mu(+) (m_mu_c2) = 105.7 MeV
    - The neutrino (nu) is assumed to be massless.

    The provided answer is D, which corresponds to:
    - KE of muon = 4.12 MeV
    - KE of neutrino = 29.8 MeV
    """

    # --- Define constants and the proposed answer ---
    m_pi_c2 = 139.6  # MeV
    m_mu_c2 = 105.7  # MeV
    
    # The answer to check is D: 4.12 MeV, 29.8 MeV
    # By convention, the lighter particle (neutrino) gets more kinetic energy.
    # However, the problem is about a two-body decay, where energies are fixed.
    # The massive particle (muon) will have the lower kinetic energy.
    ke_mu_ans = 4.12  # MeV
    ke_nu_ans = 29.8  # MeV

    # --- Verification Step 1: Direct Calculation ---
    # We can calculate the exact kinetic energies from first principles.
    # From conservation of energy and momentum, we can derive the formulas:
    # KE_nu = (m_pi_c2^2 - m_mu_c2^2) / (2 * m_pi_c2)
    # KE_mu = (m_pi_c2 - m_mu_c2) - KE_nu

    try:
        # Calculate theoretical KE of the neutrino
        ke_nu_calc = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
        
        # Calculate the total kinetic energy released (Q-value)
        q_value = m_pi_c2 - m_mu_c2
        
        # Calculate theoretical KE of the muon
        ke_mu_calc = q_value - ke_nu_calc
    except Exception as e:
        return f"An error occurred during calculation: {e}"

    # Compare the calculated values with the answer's values.
    # A tolerance is used for floating-point comparison, as the answer is rounded.
    tolerance = 0.01  # The answer is given to two decimal places.

    if not math.isclose(ke_mu_ans, ke_mu_calc, abs_tol=tolerance):
        return (f"Incorrect. The kinetic energy of the muon is wrong. "
                f"The calculated value is approximately {ke_mu_calc:.3f} MeV, "
                f"but the answer provides {ke_mu_ans} MeV.")

    if not math.isclose(ke_nu_ans, ke_nu_calc, abs_tol=tolerance):
        return (f"Incorrect. The kinetic energy of the neutrino is wrong. "
                f"The calculated value is approximately {ke_nu_calc:.3f} MeV, "
                f"but the answer provides {ke_nu_ans} MeV.")

    # --- Verification Step 2: Check Conservation Laws with Answer's Values ---
    # This confirms the answer is physically consistent.

    # Constraint 1: Conservation of Energy (Q-value)
    # The sum of the kinetic energies must equal the difference in rest masses.
    total_ke_ans = ke_mu_ans + ke_nu_ans
    if not math.isclose(total_ke_ans, q_value, abs_tol=tolerance):
        return (f"Incorrect. The answer violates conservation of energy. "
                f"The sum of the kinetic energies in the answer ({total_ke_ans:.2f} MeV) "
                f"does not match the total kinetic energy released (Q-value), which should be "
                f"{q_value:.2f} MeV.")

    # Constraint 2: Conservation of Momentum
    # The magnitudes of the momenta of the two particles must be equal.
    # For the massless neutrino: p_nu*c = KE_nu
    # For the massive muon: (p_mu*c)^2 = E_mu^2 - (m_mu*c^2)^2, where E_mu = KE_mu + m_mu*c^2
    
    p_nu_c_ans = ke_nu_ans
    
    e_mu_ans = ke_mu_ans + m_mu_c2
    p_mu_c_sq_ans = e_mu_ans**2 - m_mu_c2**2
    p_mu_c_ans = math.sqrt(p_mu_c_sq_ans)

    if not math.isclose(p_nu_c_ans, p_mu_c_ans, abs_tol=tolerance):
        return (f"Incorrect. The answer violates conservation of momentum. "
                f"The momentum of the neutrino (p*c) is {p_nu_c_ans:.2f} MeV, "
                f"but the momentum of the muon (p*c) is calculated to be {p_mu_c_ans:.2f} MeV from its KE.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)