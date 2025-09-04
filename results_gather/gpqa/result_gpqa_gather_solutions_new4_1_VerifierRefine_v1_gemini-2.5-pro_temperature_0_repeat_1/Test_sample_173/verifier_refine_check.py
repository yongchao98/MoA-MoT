import math

def check_fission_problem():
    """
    This function verifies the solution to the fission problem by recalculating all values from scratch.
    It checks the following steps:
    1. Calculation of fragment rest-mass energies.
    2. Calculation of the total kinetic energy released (Q-value).
    3. Calculation of the classical kinetic energy of the massive fragment (T1_classical).
    4. Calculation of the relativistic kinetic energy of the massive fragment (T1_relativistic).
    5. The final difference between the relativistic and classical values.
    """
    
    # Step 1: Define initial values and calculate fragment rest-mass energies
    try:
        Mc2 = 300.0  # Initial rest-mass energy in GeV
        
        # From the problem statement:
        # m1 + m2 = 0.99 * M  => m1c^2 + m2c^2 = 0.99 * Mc^2
        # m1 = 2 * m2         => m1c^2 = 2 * m2c^2
        
        # Solving the system of equations:
        # 2*m2c^2 + m2c^2 = 0.99 * Mc^2
        # 3*m2c^2 = 0.99 * Mc^2
        m2c2 = (0.99 * Mc2) / 3.0  # Rest-mass energy of the lighter fragment
        m1c2 = 2.0 * m2c2          # Rest-mass energy of the more massive fragment
        
        # Expected values: m2c2 = 99 GeV, m1c2 = 198 GeV
        if not (math.isclose(m2c2, 99.0) and math.isclose(m1c2, 198.0)):
            return f"Constraint check failed: Incorrect fragment rest-mass energies. Calculated m1c2={m1c2:.2f} GeV, m2c2={m2c2:.2f} GeV. Expected 198 GeV and 99 GeV."

    except Exception as e:
        return f"An error occurred during Step 1 (calculating fragment energies): {e}"

    # Step 2: Calculate the total kinetic energy released (Q-value)
    try:
        Q = Mc2 - (m1c2 + m2c2)
        
        # Expected value: Q = 300 - (198 + 99) = 3 GeV
        if not math.isclose(Q, 3.0):
            return f"Constraint check failed: Incorrect Q-value. Calculated Q={Q:.2f} GeV. Expected 3 GeV."
    except Exception as e:
        return f"An error occurred during Step 2 (calculating Q-value): {e}"

    # Step 3: Calculate T1 using the classical (non-relativistic) approximation
    try:
        # From conservation of momentum (p1=p2) and T = p^2/(2m), we get T1/T2 = m2/m1 = 1/2.
        # So, T2_cl = 2 * T1_cl.
        # T1_cl + T2_cl = Q => T1_cl + 2*T1_cl = Q => 3*T1_cl = Q
        T1_classical = Q / 3.0
        
        # Expected value: T1_classical = 3 / 3 = 1 GeV
        if not math.isclose(T1_classical, 1.0):
            return f"Calculation error: Incorrect classical T1. Calculated T1_classical={T1_classical:.4f} GeV. Expected 1.0 GeV."
    except Exception as e:
        return f"An error occurred during Step 3 (calculating classical T1): {e}"

    # Step 4: Calculate the correct (relativistic) T1
    try:
        # From conservation of momentum (p1=p2) and (pc)^2 = T^2 + 2T(mc^2), we get:
        # T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
        # Substitute T2 = Q - T1 and solve for T1. The simplified solution is:
        # T1_rel = (Q * (Q + 2*m2c2)) / (2 * Mc2)
        T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)
        
        # Expected value: T1_rel = (3 * (3 + 2*99)) / (2*300) = (3 * 201) / 600 = 603/600 = 1.005 GeV
        if not math.isclose(T1_relativistic, 1.005):
            return f"Calculation error: Incorrect relativistic T1. Calculated T1_relativistic={T1_relativistic:.4f} GeV. Expected 1.005 GeV."
    except Exception as e:
        return f"An error occurred during Step 4 (calculating relativistic T1): {e}"

    # Step 5: Find the difference and check against the answer
    try:
        difference_GeV = T1_relativistic - T1_classical
        difference_MeV = difference_GeV * 1000.0
        
        # The expected difference is 5 MeV, which corresponds to option A.
        expected_difference_MeV = 5.0
        
        if math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
            # The calculation confirms the result of 5 MeV.
            # The provided answer also concludes 5 MeV and selects option A.
            return "Correct"
        else:
            return f"The final calculated difference is {difference_MeV:.4f} MeV, which does not match the expected answer of 5 MeV (Option A)."
    except Exception as e:
        return f"An error occurred during Step 5 (calculating the final difference): {e}"

# Run the check
result = check_fission_problem()
print(result)