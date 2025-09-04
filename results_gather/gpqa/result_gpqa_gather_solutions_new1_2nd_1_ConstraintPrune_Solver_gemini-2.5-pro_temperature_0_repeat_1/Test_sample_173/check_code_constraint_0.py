import math

def check_correctness():
    """
    This function verifies the step-by-step calculation for the given physics problem.
    It checks if the final answer 'C' (5 MeV) is correct.
    """
    # --- Given values and constraints ---
    E_initial_GeV = 300.0
    m_sum_ratio = 0.99
    m_ratio = 2.0  # m1 / m2
    options = {'A': 10.0, 'B': 20.0, 'C': 5.0, 'D': 2.0}
    final_answer_option = 'C'

    # --- Step 1: Calculate Fragment Rest-Mass Energies ---
    # We have the system of equations:
    # 1) m1c^2 + m2c^2 = 0.99 * E_initial
    # 2) m1c^2 = 2 * m2c^2
    # Substituting (2) into (1) gives: 3 * m2c^2 = 0.99 * E_initial
    m2c2_GeV = m_sum_ratio * E_initial_GeV / (m_ratio + 1)
    m1c2_GeV = m_ratio * m2c2_GeV

    expected_m2c2 = 99.0
    expected_m1c2 = 198.0
    if not (math.isclose(m2c2_GeV, expected_m2c2) and math.isclose(m1c2_GeV, expected_m1c2)):
        return f"Constraint failed: Fragment rest-mass energies are incorrect. Calculated m1c^2={m1c2_GeV:.2f} GeV, m2c^2={m2c2_GeV:.2f} GeV."

    # --- Step 2: Calculate Total Kinetic Energy (Q-value) ---
    Q_GeV = E_initial_GeV - (m1c2_GeV + m2c2_GeV)
    expected_Q = 3.0
    if not math.isclose(Q_GeV, expected_Q):
        return f"Constraint failed: Q-value is incorrect. Calculated Q={Q_GeV:.2f} GeV."

    # --- Step 3: Calculate Classical T1 ---
    # From conservation of momentum, T1/T2 = m2/m1 = 1/2.
    # T1 + T2 = Q => T1 + 2*T1 = Q => 3*T1 = Q
    T1_classical_GeV = Q_GeV / (m_ratio + 1)
    expected_T1_classical = 1.0
    if not math.isclose(T1_classical_GeV, expected_T1_classical):
        return f"Constraint failed: Classical T1 is incorrect. Calculated T1_classical={T1_classical_GeV:.4f} GeV."

    # --- Step 4: Calculate Relativistic T1 ---
    # The relativistic calculation leads to the equation: 600 * T1 = 603
    T1_relativistic_GeV = 603.0 / 600.0
    expected_T1_relativistic = 1.005
    if not math.isclose(T1_relativistic_GeV, expected_T1_relativistic):
        return f"Constraint failed: Relativistic T1 is incorrect. Calculated T1_relativistic={T1_relativistic_GeV:.4f} GeV."

    # --- Step 5: Calculate the Difference in MeV ---
    difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    difference_MeV = difference_GeV * 1000.0
    expected_difference_MeV = 5.0
    if not math.isclose(difference_MeV, expected_difference_MeV):
        return f"Constraint failed: Final difference is incorrect. Calculated difference={difference_MeV:.2f} MeV."

    # --- Step 6: Check if the final answer option matches the calculated value ---
    if final_answer_option not in options:
        return f"Error: The provided answer option '{final_answer_option}' is not a valid choice."

    final_answer_value = options[final_answer_option]
    if not math.isclose(final_answer_value, difference_MeV):
        return f"The final answer is incorrect. The calculated difference is {difference_MeV:.2f} MeV, but the chosen option '{final_answer_option}' corresponds to {final_answer_value} MeV."

    return "Correct"

# Execute the check and print the result
result = check_correctness()
print(result)