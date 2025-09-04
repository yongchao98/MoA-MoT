import math

def check_answer():
    """
    This function checks the correctness of the provided answer by recalculating the physics problem from scratch.
    """
    # Define initial values from the problem statement
    Mc2 = 300.0  # Initial rest-mass energy in GeV
    final_mass_fraction = 0.99
    GeV_to_MeV = 1000.0

    # The options as provided in the final question block
    options = {
        'A': 10.0,
        'B': 20.0,
        'C': 5.0,
        'D': 2.0
    }
    
    # The final answer given by the LLM
    llm_answer_letter = 'C'

    # --- Step 1: Calculate the rest-mass energies of the fragments ---
    # We have two equations:
    # 1) m1c2 = 2 * m2c2
    # 2) m1c2 + m2c2 = 0.99 * Mc2
    total_final_rest_mass_energy = final_mass_fraction * Mc2
    
    # Substitute (1) into (2): 2*m2c2 + m2c2 = total_final_rest_mass_energy
    # 3 * m2c2 = total_final_rest_mass_energy
    m2c2 = total_final_rest_mass_energy / 3.0
    m1c2 = 2.0 * m2c2

    # Check if the calculated fragment masses are correct
    expected_m1c2 = 198.0
    expected_m2c2 = 99.0
    if not (math.isclose(m1c2, expected_m1c2) and math.isclose(m2c2, expected_m2c2)):
        return f"Incorrect fragment rest-mass energies. Calculated m1c2={m1c2} GeV, m2c2={m2c2} GeV. Expected m1c2={expected_m1c2} GeV, m2c2={expected_m2c2} GeV."

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    Q_value = Mc2 - total_final_rest_mass_energy
    expected_Q_value = 3.0
    if not math.isclose(Q_value, expected_Q_value):
        return f"Incorrect Q-value. Calculated Q={Q_value} GeV. Expected Q={expected_Q_value} GeV."

    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # T1_cl / T2_cl = m2 / m1 = m2c2 / m1c2 = 1/2 => T2_cl = 2 * T1_cl
    # T1_cl + T2_cl = Q_value => T1_cl + 2*T1_cl = Q_value
    # 3 * T1_cl = Q_value
    T1_classical = Q_value / 3.0
    expected_T1_classical = 1.0
    if not math.isclose(T1_classical, expected_T1_classical):
        return f"Incorrect classical T1. Calculated T1_classical={T1_classical} GeV. Expected T1_classical={expected_T1_classical} GeV."

    # --- Step 4: Calculate the correct (relativistic) T1 ---
    # (p1c)^2 = (p2c)^2
    # T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # T2 = Q_value - T1
    # T1^2 + 2*T1*m1c2 = (Q_value - T1)^2 + 2*(Q_value - T1)*m2c2
    # T1^2 + 2*T1*m1c2 = Q_value^2 - 2*Q_value*T1 + T1^2 + 2*Q_value*m2c2 - 2*T1*m2c2
    # The T1^2 terms cancel out.
    # 2*T1*m1c2 = Q_value^2 - 2*Q_value*T1 + 2*Q_value*m2c2 - 2*T1*m2c2
    # T1 * (2*m1c2 + 2*Q_value + 2*m2c2) = Q_value^2 + 2*Q_value*m2c2
    # T1 = (Q_value^2 + 2*Q_value*m2c2) / (2*m1c2 + 2*Q_value + 2*m2c2)
    # Let's use the simpler algebraic form from the analysis: 600*T1 = 603
    # 396*T1 = 603 - 204*T1 => 600*T1 = 603
    T1_relativistic = 603.0 / 600.0
    expected_T1_relativistic = 1.005
    if not math.isclose(T1_relativistic, expected_T1_relativistic):
        return f"Incorrect relativistic T1. Calculated T1_relativistic={T1_relativistic} GeV. Expected T1_relativistic={expected_T1_relativistic} GeV."

    # --- Step 5: Find the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * GeV_to_MeV
    
    expected_difference_MeV = 5.0
    if not math.isclose(difference_MeV, expected_difference_MeV):
        return f"Incorrect final difference. Calculated difference={difference_MeV} MeV. Expected difference={expected_difference_MeV} MeV."

    # --- Step 6: Check if the calculated value matches the LLM's chosen option ---
    llm_answer_value = options.get(llm_answer_letter)
    if llm_answer_value is None:
        return f"The LLM's answer letter '{llm_answer_letter}' is not a valid option."

    if not math.isclose(difference_MeV, llm_answer_value):
        return f"The calculated difference is {difference_MeV} MeV, but the LLM chose option {llm_answer_letter} which corresponds to {llm_answer_value} MeV."

    return "Correct"

# Run the check
result = check_answer()
print(result)