import math

def check_fission_energy_difference():
    """
    Checks the correctness of the calculated difference between relativistic and classical
    kinetic energy for a fission fragment.

    The problem states:
    - Initial rest-mass energy (Mc^2) = 300 GeV.
    - Final fragments m1 and m2, with m1 = 2 * m2.
    - Sum of final rest-masses (m1 + m2) = 0.99 * M.
    - T1 is the kinetic energy of the more massive fragment (m1).
    - The goal is to find the difference between the correct (relativistic) T1 and the
      classical T1.
    - The expected answer is 5 MeV.
    """
    # Initial values
    Mc2 = 300.0  # GeV

    # Step 1: Calculate the rest-mass energies of the fragments
    # We have a system of two equations:
    # 1) m1c2 + m2c2 = 0.99 * Mc2
    # 2) m1c2 = 2 * m2c2
    # Substitute (2) into (1): 2*m2c2 + m2c2 = 0.99 * Mc2
    # 3*m2c2 = 0.99 * Mc2
    m2c2 = 0.99 * Mc2 / 3.0
    m1c2 = 2.0 * m2c2

    # Check if the sum is correct
    if not math.isclose(m1c2 + m2c2, 0.99 * Mc2):
        return f"Constraint check failed: The sum of fragment rest-mass energies ({m1c2 + m2c2} GeV) does not equal 99% of the initial energy ({0.99 * Mc2} GeV)."

    # Step 2: Calculate the total kinetic energy released (Q-value)
    Q = Mc2 - (m1c2 + m2c2)

    # Step 3: Calculate T1 using the classical (non-relativistic) approximation
    # From conservation of momentum (p1=p2) and T = p^2/(2m), we get T1/T2 = m2/m1 = 1/2.
    # So, T2_classical = 2 * T1_classical.
    # T1_classical + T2_classical = Q  =>  T1_classical + 2*T1_classical = Q
    # 3 * T1_classical = Q
    T1_classical = Q / 3.0

    # Step 4: Calculate the correct (relativistic) T1 value
    # From conservation of momentum (p1=p2) and (pc)^2 = T^2 + 2*T*(mc^2), we get:
    # T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
    # We know T2 = Q - T1. Substitute this into the equation.
    # T1^2 + 2*T1*m1c2 = (Q - T1)^2 + 2*(Q - T1)*m2c2
    # T1^2 + 2*T1*m1c2 = Q^2 - 2*Q*T1 + T1^2 + 2*Q*m2c2 - 2*T1*m2c2
    # Canceling T1^2 and rearranging for T1:
    # T1 * (2*m1c2 + 2*Q + 2*m2c2) = Q^2 + 2*Q*m2c2
    # T1 * (2*(m1c2 + m2c2 + Q)) = Q * (Q + 2*m2c2)
    # Since m1c2 + m2c2 + Q = Mc2, this simplifies to:
    # T1 * (2 * Mc2) = Q * (Q + 2*m2c2)
    T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)

    # Step 5: Find the difference and convert to MeV
    difference_gev = T1_relativistic - T1_classical
    difference_mev = difference_gev * 1000.0

    # The expected answer from the analysis is 5 MeV.
    expected_difference_mev = 5.0

    # Check if the calculated difference matches the expected answer
    if math.isclose(difference_mev, expected_difference_mev, rel_tol=1e-6):
        return "Correct"
    else:
        # Provide a detailed reason for the incorrectness
        reason = (
            f"The final calculated difference is incorrect.\n"
            f"Calculated values:\n"
            f"  - m1c^2 = {m1c2:.3f} GeV\n"
            f"  - m2c^2 = {m2c2:.3f} GeV\n"
            f"  - Q-value = {Q:.3f} GeV\n"
            f"  - T1 (classical) = {T1_classical:.4f} GeV\n"
            f"  - T1 (relativistic) = {T1_relativistic:.4f} GeV\n"
            f"  - Calculated Difference = {difference_mev:.4f} MeV\n"
            f"  - Expected Difference = {expected_difference_mev} MeV\n"
            f"The calculated difference does not match the expected value of 5 MeV."
        )
        return reason

# Run the check
result = check_fission_energy_difference()
print(result)