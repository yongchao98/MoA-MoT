import sympy
import math

def check_answer():
    """
    This function checks the correctness of the provided LLM answer for the physics problem.
    It recalculates all values from scratch and compares them to the steps in the answer.
    """
    try:
        # Define constants and initial values from the problem statement
        Mc2 = 300.0  # GeV
        sum_m_ratio = 0.99
        m1_div_m2_ratio = 2.0
        GeV_to_MeV = 1000.0

        # --- Step 1: Calculate Fragment Rest-Mass Energies ---
        # m1c2 + m2c2 = 0.99 * Mc2
        # m1c2 = 2 * m2c2
        # => 3 * m2c2 = 0.99 * Mc2
        m1c2_plus_m2c2 = sum_m_ratio * Mc2
        
        # 2*m2c2 + m2c2 = m1c2_plus_m2c2
        # 3*m2c2 = m1c2_plus_m2c2
        m2c2 = m1c2_plus_m2c2 / (m1_div_m2_ratio + 1)
        m1c2 = m1_div_m2_ratio * m2c2

        # Verification for Step 1
        expected_m1c2 = 198.0
        expected_m2c2 = 99.0
        if not (math.isclose(m1c2, expected_m1c2) and math.isclose(m2c2, expected_m2c2)):
            return f"Incorrect fragment rest-mass energies. Calculated m1c^2={m1c2} GeV, m2c^2={m2c2} GeV. Expected m1c^2={expected_m1c2} GeV, m2c^2={expected_m2c2} GeV."

        # --- Step 2: Calculate Total Kinetic Energy Released (Q-value) ---
        Q = Mc2 - m1c2_plus_m2c2
        
        # Verification for Step 2
        expected_Q = 3.0
        if not math.isclose(Q, expected_Q):
            return f"Incorrect Q-value. Calculated Q={Q} GeV. Expected Q={expected_Q} GeV."

        # --- Step 3: Calculate T1 using the Classical (Non-Relativistic) Approximation ---
        # T1_cl / T2_cl = m2 / m1 = 1/2 => T2_cl = 2 * T1_cl
        # T1_cl + T2_cl = Q
        # T1_cl + 2*T1_cl = Q => 3*T1_cl = Q
        T1_classical = Q / 3.0
        
        # Verification for Step 3
        expected_T1_classical = 1.0
        if not math.isclose(T1_classical, expected_T1_classical):
            return f"Incorrect classical T1. Calculated T1_classical={T1_classical} GeV. Expected T1_classical={expected_T1_classical} GeV."

        # --- Step 4: Calculate the Correct (Relativistic) T1 ---
        # Using sympy to solve the equation to avoid manual algebra errors
        # (p1c)^2 = (p2c)^2
        # T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
        # T2 = Q - T1
        T1_sym = sympy.Symbol('T1_sym')
        T2_sym = Q - T1_sym
        
        # Relativistic momentum-energy relation: (pc)^2 = T^2 + 2*T*(mc^2)
        eq = sympy.Eq(T1_sym**2 + 2 * T1_sym * m1c2, T2_sym**2 + 2 * T2_sym * m2c2)
        
        # Solve for T1
        solutions = sympy.solve(eq, T1_sym)
        # The solution should be a single positive value
        if len(solutions) != 1:
            return f"Relativistic equation solving failed. Found {len(solutions)} solutions: {solutions}."
        T1_relativistic = float(solutions[0])

        # Verification for Step 4
        expected_T1_relativistic = 1.005
        if not math.isclose(T1_relativistic, expected_T1_relativistic, rel_tol=1e-6):
            return f"Incorrect relativistic T1. Calculated T1_relativistic={T1_relativistic} GeV. Expected T1_relativistic={expected_T1_relativistic} GeV."

        # --- Step 5: Find the Difference and Convert Units ---
        difference_GeV = T1_relativistic - T1_classical
        difference_MeV = difference_GeV * GeV_to_MeV

        # Verification for Step 5
        expected_difference_MeV = 5.0
        if not math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
            return f"Incorrect final difference. Calculated difference={difference_MeV} MeV. Expected difference={expected_difference_MeV} MeV."

        # --- Final Check: Match with options ---
        # The provided answer is <<<C>>>.
        # The options are: A) 20 MeV, B) 2 MeV, C) 5 MeV, D) 10 MeV.
        # Option C corresponds to 5 MeV.
        final_answer_choice = 'C'
        options = {'A': 20, 'B': 2, 'C': 5, 'D': 10}
        
        if math.isclose(options[final_answer_choice], difference_MeV, rel_tol=1e-6):
            return "Correct"
        else:
            return f"The calculated difference is {difference_MeV} MeV, which corresponds to option C. However, the provided answer choice '{final_answer_choice}' does not match the calculated value in the context of the options."

    except Exception as e:
        return f"An error occurred during the check: {e}"

# Run the check
result = check_answer()
print(result)