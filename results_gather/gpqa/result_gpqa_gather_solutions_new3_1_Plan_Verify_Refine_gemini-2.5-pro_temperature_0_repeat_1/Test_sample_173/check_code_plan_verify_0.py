import math

def check_fission_energy_difference():
    """
    This function checks the correctness of the answer to the physics problem.
    It calculates the difference between the relativistic and classical kinetic energy
    of the more massive fragment and compares it to the provided answer.
    """
    # 1. Define initial conditions and constants
    E_M = 300.0  # Initial rest-mass energy in GeV

    # 2. Calculate the rest-mass energies of the fragments
    # m1c^2 + m2c^2 = 0.99 * E_M
    # m1c^2 = 2 * m2c^2
    m1c2_plus_m2c2 = 0.99 * E_M
    # 2*m2c^2 + m2c^2 = m1c2_plus_m2c2 => 3*m2c^2 = m1c2_plus_m2c2
    m2c2 = m1c2_plus_m2c2 / 3.0
    m1c2 = 2.0 * m2c2

    # 3. Calculate the total kinetic energy released (Q-value)
    # Q = E_M - (m1c^2 + m2c^2)
    Q = E_M - m1c2_plus_m2c2
    
    # Check if the calculated values match the reasoning in the provided answers
    if not math.isclose(m1c2, 198.0) or not math.isclose(m2c2, 99.0) or not math.isclose(Q, 3.0):
        return f"Intermediate calculation mismatch. Calculated m1c^2={m1c2}, m2c^2={m2c2}, Q={Q}. Expected 198, 99, 3 respectively."

    # 4. Calculate T1 using the classical (non-relativistic) approximation
    # From conservation of momentum, T1_cl / T2_cl = m2 / m1 = 1/2 => T2_cl = 2 * T1_cl
    # T1_cl + T2_cl = Q => T1_cl + 2*T1_cl = Q => 3*T1_cl = Q
    T1_classical = Q / 3.0

    # 5. Calculate T1 using the correct (relativistic) mechanics
    # From conservation of momentum, (p1c)^2 = (p2c)^2
    # (pc)^2 = T^2 + 2*T*m*c^2
    # T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # Substitute T2 = Q - T1 and solve the resulting linear equation for T1:
    # T1^2 + 2*T1*m1c2 = (Q - T1)^2 + 2*(Q - T1)*m2c2
    # T1^2 + 2*T1*m1c2 = Q^2 - 2*Q*T1 + T1^2 + 2*Q*m2c2 - 2*T1*m2c2
    # 2*T1*m1c2 = Q^2 - 2*Q*T1 + 2*Q*m2c2 - 2*T1*m2c2
    # T1 * (2*m1c2 + 2*Q + 2*m2c2) = Q^2 + 2*Q*m2c2
    # T1 * (2 * (m1c2 + m2c2 + Q)) = Q * (Q + 2*m2c2)
    # T1 * (2 * E_M) = Q * (Q + 2*m2c2)
    T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * E_M)

    # Check if the calculated T1 values match the reasoning
    if not math.isclose(T1_classical, 1.0):
        return f"Classical T1 calculation is incorrect. Calculated T1_classical = {T1_classical} GeV, expected 1.0 GeV."
    if not math.isclose(T1_relativistic, 1.005):
        return f"Relativistic T1 calculation is incorrect. Calculated T1_relativistic = {T1_relativistic} GeV, expected 1.005 GeV."

    # 6. Find the difference and convert to MeV
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * 1000.0

    # 7. Check against the provided answer
    # The provided answer is <<<B>>>, which corresponds to 5 MeV.
    expected_value_MeV = 5.0
    
    if math.isclose(difference_MeV, expected_value_MeV, rel_tol=1e-5):
        return "Correct"
    else:
        return (f"The calculated difference is {difference_MeV:.4f} MeV, "
                f"which does not match the expected value of {expected_value_MeV} MeV from option B.")

# Run the check
result = check_fission_energy_difference()
print(result)