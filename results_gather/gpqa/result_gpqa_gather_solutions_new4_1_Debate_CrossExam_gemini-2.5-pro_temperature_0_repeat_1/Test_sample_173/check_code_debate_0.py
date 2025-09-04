import math

def check_answer():
    """
    Checks the correctness of the provided answer to the physics problem.
    """
    # Define constants and initial values from the problem statement
    # All energy values are in GeV unless otherwise specified.
    Mc2 = 300.0  # Rest-mass energy of the initial nucleus

    # The final answer from the LLM is 'B', which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # --- Step 1: Calculate the rest-mass energies of the fragments ---
    # We have two equations:
    # 1) m1 = 2 * m2  => m1c2 = 2 * m2c2
    # 2) m1 + m2 = 0.99 * M => m1c2 + m2c2 = 0.99 * Mc2
    sum_m_c2 = 0.99 * Mc2
    
    # Substitute (1) into (2): 2*m2c2 + m2c2 = sum_m_c2
    # 3 * m2c2 = sum_m_c2
    m2c2 = sum_m_c2 / 3.0
    m1c2 = 2.0 * m2c2

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # Q = Initial rest energy - Final rest energy
    Q_value = Mc2 - (m1c2 + m2c2)

    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # In classical mechanics, T = p^2 / (2m).
    # Since momentum p is conserved (p1=p2), T1/T2 = m2/m1 = 1/2.
    # So, T2_classical = 2 * T1_classical.
    # T1_classical + T2_classical = Q_value
    # T1_classical + 2 * T1_classical = Q_value
    # 3 * T1_classical = Q_value
    T1_classical = Q_value / 3.0

    # --- Step 4: Calculate the correct (relativistic) T1 ---
    # The relativistic relation is (pc)^2 = T^2 + 2*T*(mc^2).
    # Since momentum p is conserved (p1=p2), we have:
    # T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # We also know T2 = Q_value - T1.
    # T1^2 + 2*T1*m1c2 = (Q_value - T1)^2 + 2*(Q_value - T1)*m2c2
    # T1^2 + 2*T1*m1c2 = Q_value^2 - 2*Q_value*T1 + T1^2 + 2*Q_value*m2c2 - 2*T1*m2c2
    # The T1^2 terms cancel. Rearranging to solve for T1:
    # 2*T1*m1c2 + 2*Q_value*T1 + 2*T1*m2c2 = Q_value^2 + 2*Q_value*m2c2
    # T1 * (2*m1c2 + 2*Q_value + 2*m2c2) = Q_value * (Q_value + 2*m2c2)
    # T1 * (2 * (m1c2 + m2c2 + Q_value)) = Q_value * (Q_value + 2*m2c2)
    # Since m1c2 + m2c2 + Q_value = Mc2, the equation simplifies to:
    # T1 * (2 * Mc2) = Q_value * (Q_value + 2*m2c2)
    T1_relativistic = (Q_value * (Q_value + 2 * m2c2)) / (2 * Mc2)

    # --- Step 5: Find the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    calculated_difference_MeV = difference_GeV * 1000.0

    # --- Step 6: Check correctness ---
    # Use a small tolerance for floating-point comparison
    tolerance = 1e-6
    if abs(calculated_difference_MeV - expected_difference_MeV) < tolerance:
        return "Correct"
    else:
        # Provide a detailed reason for the error
        reason = (
            f"The final answer is incorrect.\n"
            f"The provided answer corresponds to a difference of {expected_difference_MeV} MeV.\n"
            f"However, the calculation shows a different result.\n\n"
            f"Here are the calculated values:\n"
            f"- Rest-mass energy of massive fragment (m1c2): {m1c2:.3f} GeV\n"
            f"- Rest-mass energy of lighter fragment (m2c2): {m2c2:.3f} GeV\n"
            f"- Total kinetic energy released (Q-value): {Q_value:.3f} GeV\n"
            f"- Classical T1: {T1_classical:.3f} GeV\n"
            f"- Relativistic T1: {T1_relativistic:.3f} GeV\n"
            f"- Calculated difference: {T1_relativistic:.3f} GeV - {T1_classical:.3f} GeV = {difference_GeV:.3f} GeV\n"
            f"- Calculated difference in MeV: {calculated_difference_MeV:.3f} MeV\n"
            f"This calculated value does not match the expected value of {expected_difference_MeV} MeV."
        )
        return reason

# Run the check
result = check_answer()
print(result)