import math

def check_correctness():
    """
    Checks the correctness of the answer for the pion decay problem.

    The problem asks for the kinetic energies (KE) of the product particles
    in the decay: Pi(+) -> mu(+) + nu, where the Pi(+) is stationary.

    Given:
    - Rest mass energy of Pi(+) (m_pi_c2): 139.6 MeV
    - Rest mass energy of mu(+) (m_mu_c2): 105.7 MeV
    - The neutrino (nu) is assumed to be massless.

    The provided answer is 'B', which corresponds to the values 4.12 MeV and 29.8 MeV.
    """

    # --- Given values and the answer to check ---
    m_pi_c2 = 139.6  # MeV
    m_mu_c2 = 105.7  # MeV
    
    # The final answer provided is 'B', with values 4.12 MeV and 29.8 MeV.
    # In a two-body decay from rest, the lighter particle (neutrino) gets more kinetic energy.
    answer_values = sorted([4.12, 29.8])
    ans_ke_mu = answer_values[0]
    ans_ke_nu = answer_values[1]

    # --- Calculation from first principles ---

    # 1. Calculate the total kinetic energy released (Q-value).
    # This must be conserved.
    q_value = m_pi_c2 - m_mu_c2
    
    # 2. Calculate the kinetic energy of the massless neutrino (KE_nu).
    # This formula is derived from conservation of energy and momentum.
    try:
        calc_ke_nu = (m_pi_c2**2 - m_mu_c2**2) / (2 * m_pi_c2)
    except ZeroDivisionError:
        return "Calculation Error: The rest mass of the pion cannot be zero."

    # 3. Calculate the kinetic energy of the massive muon (KE_mu).
    calc_ke_mu = q_value - calc_ke_nu

    # --- Verification ---
    
    # Set a tolerance for comparing floating-point numbers.
    # A relative tolerance of 1% is reasonable given the rounding in the options.
    tolerance = 0.01

    # Constraint 1: The sum of the kinetic energies in the answer must equal the Q-value.
    if not math.isclose(ans_ke_mu + ans_ke_nu, q_value, rel_tol=tolerance):
        return (f"Incorrect. The sum of the kinetic energies in the answer ({ans_ke_mu + ans_ke_nu:.2f} MeV) "
                f"does not match the total kinetic energy released (Q-value), which should be "
                f"{q_value:.2f} MeV. This violates the conservation of energy.")

    # Constraint 2: The individual calculated kinetic energies must match the answer's values.
    if not math.isclose(calc_ke_mu, ans_ke_mu, rel_tol=tolerance):
        return (f"Incorrect. The calculated kinetic energy for the muon is {calc_ke_mu:.2f} MeV, "
                f"but the answer provides {ans_ke_mu:.2f} MeV.")

    if not math.isclose(calc_ke_nu, ans_ke_nu, rel_tol=tolerance):
        return (f"Incorrect. The calculated kinetic energy for the neutrino is {calc_ke_nu:.2f} MeV, "
                f"but the answer provides {ans_ke_nu:.2f} MeV.")

    # If all checks pass, the answer is correct.
    return "Correct"

# Execute the check
result = check_correctness()
print(result)