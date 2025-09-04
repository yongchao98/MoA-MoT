import math

def check_correctness_of_answer():
    """
    This function checks the correctness of the provided LLM answer by recalculating the physics problem from scratch.
    The LLM's final answer is 'D', which corresponds to 5 MeV. This function will verify if that value is correct.
    """
    
    # --- Problem Setup ---
    # The final answer provided by the LLM is 'D'.
    # The options are: A) 2 MeV, B) 20 MeV, C) 10 MeV, D) 5 MeV.
    # Therefore, the value to check against is 5 MeV.
    expected_difference_MeV = 5.0
    
    # --- Calculation from First Principles ---
    
    # Given values
    E_M = 300.0  # Initial rest-mass energy in GeV

    # Step 1: Calculate rest-mass energies of the fragments
    # From the problem statement:
    # m1 = 2 * m2
    # m1 + m2 = 0.99 * M
    # Substituting the first into the second: 2*m2 + m2 = 0.99*M => 3*m2 = 0.99*M => m2 = 0.33*M
    # Therefore, m1 = 2 * 0.33*M = 0.66*M.
    # We can apply these ratios to the rest-mass energies.
    m1c2 = 0.66 * E_M  # Rest-mass energy of the more massive fragment (GeV)
    m2c2 = 0.33 * E_M  # Rest-mass energy of the less massive fragment (GeV)

    # Step 2: Calculate the total kinetic energy released (Q-value)
    # The Q-value is the energy released from the conversion of rest mass.
    # Q = E_M - (m1c2 + m2c2)
    Q = E_M - (m1c2 + m2c2)

    # Step 3: Calculate T1 using the classical (non-relativistic) approximation
    # From conservation of momentum (p1=p2), the kinetic energies are inversely proportional to mass: T1_cl / T2_cl = m2 / m1 = 1/2.
    # With the total kinetic energy being T1_cl + T2_cl = Q, we get T1_cl + 2*T1_cl = Q => 3*T1_cl = Q.
    T1_classical_GeV = Q / 3.0

    # Step 4: Calculate T1 using the correct (relativistic) formula
    # From conservation of momentum, (p1c)^2 = (p2c)^2.
    # Using the relativistic relation (pc)^2 = T^2 + 2*T*(mc^2) and T2 = Q - T1, we get:
    # T1^2 + 2*T1*m1c2 = (Q - T1)^2 + 2*(Q - T1)*m2c2
    # Expanding and simplifying this equation leads to:
    # T1 * (2 * E_M) = Q * (Q + 2*m2c2)
    T1_relativistic_GeV = (Q * (Q + 2 * m2c2)) / (2 * E_M)

    # Step 5: Calculate the difference and convert to MeV
    calculated_difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    calculated_difference_MeV = calculated_difference_GeV * 1000.0

    # --- Verification ---
    # Check if the calculated value matches the value from the LLM's answer.
    # A small tolerance is used for floating-point comparison.
    if math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return "Correct"
    else:
        # If the calculation does not match, provide a detailed reason.
        reason = (
            f"Incorrect. The provided answer 'D' corresponds to a difference of {expected_difference_MeV} MeV. "
            f"However, the calculation from first principles yields a different result.\n"
            f"Calculated values:\n"
            f"  - Rest-mass energy of massive fragment (m1c2): {m1c2:.2f} GeV\n"
            f"  - Rest-mass energy of light fragment (m2c2): {m2c2:.2f} GeV\n"
            f"  - Total kinetic energy (Q): {Q:.2f} GeV\n"
            f"  - Classical T1: {T1_classical_GeV:.4f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic_GeV:.4f} GeV\n"
            f"  - Calculated difference: {calculated_difference_MeV:.4f} MeV."
        )
        return reason

# Execute the check
result = check_correctness_of_answer()
print(result)