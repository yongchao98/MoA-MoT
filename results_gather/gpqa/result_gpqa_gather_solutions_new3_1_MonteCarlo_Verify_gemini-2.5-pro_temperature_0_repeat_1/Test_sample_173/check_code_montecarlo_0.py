import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It recalculates all values from the problem statement and compares the final result
    to the given answer.

    The problem asks for the difference between the relativistic and classical kinetic
    energy of the more massive fragment.
    """
    
    # --- Step 1: Define initial values from the problem statement ---
    # The initial rest-mass energy of the nucleus is given in GeV.
    # We use floating-point numbers for precision.
    try:
        Mc2 = 300.0  # GeV
    except Exception as e:
        return f"Failed to initialize constants: {e}"

    # --- Step 2: Calculate the rest-mass energies of the two fragments ---
    # The sum of the fragments' rest-mass energies is 99% of the initial energy.
    # m1c2 + m2c2 = 0.99 * Mc2
    # The more massive fragment (m1) is twice the mass of the lighter one (m2).
    # m1c2 = 2 * m2c2
    # We can solve this system of two linear equations:
    # Substitute the second equation into the first: 2*m2c2 + m2c2 = 0.99 * Mc2
    # This gives: 3 * m2c2 = 0.99 * Mc2
    try:
        m2c2 = (0.99 * Mc2) / 3.0
        m1c2 = 2.0 * m2c2
    except Exception as e:
        return f"An error occurred during the calculation of fragment masses: {e}"

    # --- Step 3: Calculate the total kinetic energy released (Q-value) ---
    # The Q-value is the difference between the initial and final rest-mass energies.
    # This energy is converted into the kinetic energy of the fragments.
    # Q = T1 + T2 = Mc2 - (m1c2 + m2c2)
    Q = Mc2 - (m1c2 + m2c2)
    
    # Sanity check: Q should be 0.01 * Mc2 = 3.0 GeV.
    if not math.isclose(Q, 3.0, rel_tol=1e-9):
        return f"Incorrect Q-value calculation. Expected 3.0 GeV, but got {Q} GeV."

    # --- Step 4: Calculate T1 using the classical (non-relativistic) approximation ---
    # In classical mechanics, from conservation of momentum (p1=p2) and T = p^2/(2m),
    # we get T1/T2 = m2/m1 = 1/2, so T2 = 2*T1.
    # Since T1 + T2 = Q, we have T1 + 2*T1 = Q, which means 3*T1 = Q.
    T1_classical = Q / 3.0
    
    # Check this intermediate value against the provided solution's steps.
    if not math.isclose(T1_classical, 1.0, rel_tol=1e-9):
        return f"Incorrect classical T1 calculation. The provided solution calculates it as 1.0 GeV, but the code got {T1_classical} GeV."

    # --- Step 5: Calculate the correct (relativistic) T1 ---
    # Using conservation of momentum (p1=p2) and the relativistic relation (pc)^2 = T^2 + 2*T*(mc^2),
    # we get: T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2).
    # Substitute T2 = Q - T1 and solve for T1.
    # A simplified final formula derived from this is: T1 = (Q * (Q + 2*m2c2)) / (2 * Mc2)
    T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)

    # Check this intermediate value against the provided solution's steps.
    if not math.isclose(T1_relativistic, 1.005, rel_tol=1e-9):
        return f"Incorrect relativistic T1 calculation. The provided solution calculates it as 1.005 GeV, but the code got {T1_relativistic} GeV."

    # --- Step 6: Find the difference and convert to MeV ---
    # The question asks for the difference between the relativistic and classical values.
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * 1000.0

    # --- Step 7: Verify the final answer ---
    # The provided answer is <<<B>>>. The options given in the question are:
    # A) 2 MeV, B) 5 MeV, C) 10 MeV, D) 20 MeV.
    # Therefore, the expected answer is 5 MeV.
    expected_difference_MeV = 5.0
    
    if math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return "Correct"
    else:
        return (f"The final calculated difference is {difference_MeV:.3f} MeV, which does not match the "
                f"expected value of {expected_difference_MeV} MeV corresponding to answer option B. "
                f"Therefore, the provided answer <<<B>>> is incorrect.")

# Run the check
result = check_correctness()
print(result)