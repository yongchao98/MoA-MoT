import math

def check_correctness():
    """
    Checks the correctness of the provided answer to the physics problem.

    The problem asks for the difference between the relativistic (correct) and classical
    kinetic energy of the more massive fragment from a nuclear fission.
    """
    
    # --- Problem Parameters ---
    # Initial rest-mass energy of the nucleus in GeV
    Mc2 = 300.0
    
    # The final answer to check, corresponding to option D (5 MeV)
    expected_difference_MeV = 5.0

    # --- Step 1: Calculate Fragment Rest-Mass Energies ---
    # The sum of the rest-mass energies of the fragments is 99% of the initial energy.
    # m1c2 + m2c2 = 0.99 * Mc2
    sum_m_c2 = 0.99 * Mc2
    
    # We are given m1 = 2 * m2, which implies m1c2 = 2 * m2c2.
    # Substitute this into the sum: 2*m2c2 + m2c2 = sum_m_c2
    # 3 * m2c2 = sum_m_c2
    m2c2 = sum_m_c2 / 3.0
    m1c2 = 2.0 * m2c2
    
    # --- Step 2: Calculate Total Kinetic Energy Released (Q-value) ---
    # The Q-value is the mass defect converted to energy.
    Q = Mc2 - sum_m_c2
    
    # --- Step 3: Calculate T1 using Classical Approximation ---
    # In classical mechanics, T = p^2 / (2m).
    # From conservation of momentum (p1 = p2), we get T1/T2 = m2/m1 = 1/2.
    # So, T2_classical = 2 * T1_classical.
    # Since T1_classical + T2_classical = Q, we have T1_classical + 2*T1_classical = Q.
    # 3 * T1_classical = Q
    T1_classical_GeV = Q / 3.0
    
    # --- Step 4: Calculate T1 using Correct Relativistic Mechanics ---
    # The relativistic energy-momentum relation is E^2 = (pc)^2 + (mc^2)^2.
    # This can be written in terms of kinetic energy T as (pc)^2 = T^2 + 2*T*(mc^2).
    # From conservation of momentum (p1 = p2), we equate (p1c)^2 and (p2c)^2:
    # T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
    # We know T2 = Q - T1. Substituting this gives a linear equation for T1:
    # T1 * (2*m1c2 + 2*Q + 2*m2c2) = Q^2 + 2*Q*m2c2
    # The term in the parenthesis simplifies: 2*(m1c2 + m2c2 + Q) = 2*Mc2
    # So, T1 * (2 * Mc2) = Q * (Q + 2*m2c2)
    T1_correct_GeV = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)
    
    # --- Step 5: Find the Difference ---
    # Calculate the difference in GeV
    difference_GeV = T1_correct_GeV - T1_classical_GeV
    
    # Convert the result to MeV (1 GeV = 1000 MeV)
    calculated_difference_MeV = difference_GeV * 1000.0
    
    # --- Step 6: Check the Correctness of the Answer ---
    # The provided answer is <<<D>>>, which corresponds to 5 MeV.
    # We check if our calculated value matches the expected value.
    
    # Check constraints and intermediate values
    if not math.isclose(m1c2, 198.0):
        return f"Incorrect intermediate calculation: m1c^2 should be 198 GeV, but was calculated as {m1c2} GeV."
    if not math.isclose(m2c2, 99.0):
        return f"Incorrect intermediate calculation: m2c^2 should be 99 GeV, but was calculated as {m2c2} GeV."
    if not math.isclose(Q, 3.0):
        return f"Incorrect intermediate calculation: Q-value should be 3 GeV, but was calculated as {Q} GeV."
    if not math.isclose(T1_classical_GeV, 1.0):
        return f"Incorrect classical calculation: T1_classical should be 1.0 GeV, but was calculated as {T1_classical_GeV} GeV."
    if not math.isclose(T1_correct_GeV, 1.005):
        return f"Incorrect relativistic calculation: T1_correct should be 1.005 GeV, but was calculated as {T1_correct_GeV} GeV."

    # Final check on the difference
    if math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return "Correct"
    else:
        return (f"The answer is incorrect. "
                f"The calculated difference is {calculated_difference_MeV:.4f} MeV, "
                f"but the expected answer is {expected_difference_MeV} MeV.")

# Run the check
result = check_correctness()
print(result)