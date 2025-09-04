import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided LLM answer by recalculating the physics problem from scratch.
    It verifies all constraints and calculations step-by-step.
    """
    
    # --- Problem Setup ---
    # The question provides the following options:
    # A) 5 MeV.
    # B) 20 MeV.
    # C) 10 MeV.
    # D) 2 MeV.
    # The LLM's final answer is <<<A>>>.
    
    options = {'A': 5.0, 'B': 20.0, 'C': 10.0, 'D': 2.0}
    llm_answer_letter = 'A'
    expected_difference_MeV = options[llm_answer_letter]

    # --- Step-by-step Calculation ---

    # 1. Define initial constants and values from the question
    Mc2 = 300.0  # Initial rest-mass energy in GeV
    GeV_to_MeV = 1000.0

    # 2. Calculate fragment rest-mass energies based on the given constraints
    # Constraint 1: m1 = 2 * m2 -> m1c2 = 2 * m2c2
    # Constraint 2: m1 + m2 = 0.99 * M -> m1c2 + m2c2 = 0.99 * Mc2
    # We have a system of two equations. Substitute Constraint 1 into Constraint 2:
    # 2*m2c2 + m2c2 = 0.99 * Mc2 -> 3*m2c2 = 0.99 * Mc2
    m2c2 = (0.99 * Mc2) / 3.0  # Rest-mass energy of the less massive fragment
    m1c2 = 2.0 * m2c2         # Rest-mass energy of the more massive fragment

    # Verify that the constraints are met by our calculated values
    if not math.isclose(m1c2 + m2c2, 0.99 * Mc2):
        return f"Constraint check failed: The sum of fragment rest-mass energies ({m1c2 + m2c2:.2f} GeV) should be 99% of the initial energy ({0.99 * Mc2:.2f} GeV)."
    if not math.isclose(m1c2, 2 * m2c2):
        return f"Constraint check failed: The more massive fragment's energy ({m1c2:.2f} GeV) should be twice the less massive one's ({m2c2:.2f} GeV)."

    # 3. Calculate total kinetic energy released (Q-value)
    # Q is the energy released from the mass defect, which is converted to kinetic energy.
    Q = Mc2 - (m1c2 + m2c2)

    # 4. Calculate T1 using the classical (non-relativistic) approximation
    # From conservation of momentum (p1=p2), T_cl = p^2/(2m), so T1/T2 = m2/m1 = 1/2.
    # This means T2_cl = 2 * T1_cl.
    # Since T1_cl + T2_cl = Q, we have T1_cl + 2*T1_cl = Q -> 3*T1_cl = Q.
    T1_classical = Q / 3.0

    # 5. Calculate the correct (relativistic) T1
    # From conservation of momentum (p1=p2), we use the relativistic relation (pc)^2 = T^2 + 2*T*(mc^2).
    # Equating the momenta: (p1c)^2 = (p2c)^2 -> T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # Substitute T2 = Q - T1.
    # After algebraic simplification, the equation becomes: T1_rel * (2 * Mc2) = Q * (Q + 2 * m2c2)
    T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)

    # 6. Calculate the final difference and convert to MeV
    difference_GeV = T1_relativistic - T1_classical
    calculated_difference_MeV = difference_GeV * GeV_to_MeV

    # --- Verification ---
    # Compare the calculated result with the value from the LLM's answer.
    # Use a small tolerance for floating-point comparison.
    if math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return "Correct"
    else:
        reason = (
            f"The provided answer 'A' corresponds to {expected_difference_MeV} MeV, but the calculated difference is {calculated_difference_MeV:.3f} MeV.\n"
            f"The answer is incorrect.\n"
            f"Calculation Breakdown:\n"
            f"  - Initial rest-mass energy (Mc^2): {Mc2:.2f} GeV\n"
            f"  - More massive fragment rest-mass energy (m1c^2): {m1c2:.2f} GeV\n"
            f"  - Less massive fragment rest-mass energy (m2c^2): {m2c2:.2f} GeV\n"
            f"  - Total kinetic energy released (Q): {Q:.2f} GeV\n"
            f"  - Classical T1: {T1_classical:.4f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic:.4f} GeV\n"
            f"  - Difference in GeV: {difference_GeV:.4f} GeV\n"
            f"  - Difference in MeV: {calculated_difference_MeV:.3f} MeV"
        )
        return reason

# Execute the check and print the result
print(check_correctness_of_answer())