import sympy

def check_correctness():
    """
    This function checks the correctness of the provided answer by recalculating the solution from first principles.
    """
    
    # --- Problem Parameters ---
    # The question states the options are: A) 5 MeV, B) 10 MeV, C) 2 MeV, D) 20 MeV.
    # The provided answer is 'A'.
    options = {'A': 5.0, 'B': 10.0, 'C': 2.0, 'D': 20.0}
    provided_answer_letter = 'A'

    # Initial rest-mass energy in GeV
    Mc2 = 300.0

    # --- Step 1: Calculate Fragment Rest-Mass Energies ---
    # We have two equations:
    # 1) m1c2 + m2c2 = 0.99 * Mc2
    # 2) m1c2 = 2 * m2c2
    # Substituting (2) into (1): 2*m2c2 + m2c2 = 0.99 * Mc2 => 3*m2c2 = 0.99 * Mc2
    
    sum_of_final_rest_mass_energies = 0.99 * Mc2
    if not abs(sum_of_final_rest_mass_energies - 297.0) < 1e-9:
        return f"Constraint check failed: Sum of final rest-mass energies should be 297 GeV, but calculated as {sum_of_final_rest_mass_energies} GeV."

    m2c2 = sum_of_final_rest_mass_energies / 3.0
    m1c2 = 2.0 * m2c2

    # Verify the calculated rest-mass energies
    if not (abs(m1c2 - 198.0) < 1e-9 and abs(m2c2 - 99.0) < 1e-9):
        return f"Calculation error: Fragment rest-mass energies are incorrect. m1c2={m1c2} GeV, m2c2={m2c2} GeV."

    # --- Step 2: Calculate Total Kinetic Energy Released (Q-value) ---
    # Q = Initial Energy - Final Rest Energy
    Q_value = Mc2 - (m1c2 + m2c2)
    
    if not abs(Q_value - 3.0) < 1e-9:
        return f"Calculation error: Total kinetic energy (Q-value) should be 3.0 GeV, but calculated as {Q_value} GeV."

    # --- Step 3: Calculate T1 using the Classical (Non-Relativistic) Approximation ---
    # From conservation of momentum (p1=p2) and T = p^2/(2m), we get T1/T2 = m2/m1 = 1/2.
    # So, T2_cl = 2 * T1_cl.
    # T1_cl + T2_cl = Q_value => T1_cl + 2*T1_cl = Q_value => 3*T1_cl = Q_value
    T1_classical = Q_value / 3.0

    if not abs(T1_classical - 1.0) < 1e-9:
        return f"Calculation error: Classical T1 should be 1.0 GeV, but calculated as {T1_classical} GeV."

    # --- Step 4: Calculate the Correct (Relativistic) T1 ---
    # From conservation of momentum (p1=p2) and (pc)^2 = T^2 + 2*T*(mc^2), we get:
    # T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
    # We also know T2 = Q_value - T1.
    # We solve for T1.
    T1_rel = sympy.symbols('T1_rel')
    T2_rel = Q_value - T1_rel
    
    # (p1c)^2 = (p2c)^2
    equation = sympy.Eq(T1_rel**2 + 2 * T1_rel * m1c2, T2_rel**2 + 2 * T2_rel * m2c2)
    
    # Solve the equation for T1_rel
    solutions = sympy.solve(equation, T1_rel)
    
    # There should be one physically meaningful (positive) solution
    if len(solutions) != 1 or solutions[0] <= 0:
        return f"Failed to find a unique positive solution for relativistic T1. Solutions found: {solutions}"
    
    T1_relativistic = float(solutions[0])

    if not abs(T1_relativistic - 1.005) < 1e-9:
        return f"Calculation error: Relativistic T1 should be 1.005 GeV, but calculated as {T1_relativistic} GeV."

    # --- Step 5: Find the Difference and Convert to MeV ---
    difference_gev = T1_relativistic - T1_classical
    difference_mev = difference_gev * 1000.0

    # --- Step 6: Check against the provided answer ---
    expected_value = 5.0
    if abs(difference_mev - expected_value) > 1e-6:
        return f"Final calculation is incorrect. The calculated difference is {difference_mev:.3f} MeV, but the expected value from the consensus of calculations is {expected_value} MeV."

    # Check if the provided answer letter corresponds to the correct value
    if provided_answer_letter not in options:
        return f"The provided answer letter '{provided_answer_letter}' is not a valid option."
        
    selected_value = options[provided_answer_letter]
    
    if abs(selected_value - difference_mev) < 1e-6:
        return "Correct"
    else:
        return f"The provided answer is '{provided_answer_letter}', which corresponds to {selected_value} MeV. However, the calculated correct value is {difference_mev:.3f} MeV. The final answer choice is inconsistent with the calculation."

# Run the check
result = check_correctness()
print(result)