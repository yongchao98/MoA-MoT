import math

def check_answer():
    """
    Checks the correctness of the physics problem solution.
    
    The problem asks for the difference between the relativistic and classical
    kinetic energy of the more massive fragment from a fission event.
    """
    
    # 1. Define initial conditions from the problem statement
    Mc2 = 300.0  # Initial rest-mass energy in GeV
    
    # 2. Calculate fragment rest-mass energies
    # The sum of the final rest-mass energies is 99% of the initial.
    final_rest_mass_energy_total = 0.99 * Mc2
    
    # We have a system of two equations:
    # eq1: m1c2 + m2c2 = final_rest_mass_energy_total
    # eq2: m1c2 = 2 * m2c2
    # Substitute eq2 into eq1:
    # 2*m2c2 + m2c2 = final_rest_mass_energy_total
    # 3*m2c2 = final_rest_mass_energy_total
    
    m2c2 = final_rest_mass_energy_total / 3.0  # Rest-mass energy of the lighter fragment
    m1c2 = 2.0 * m2c2                         # Rest-mass energy of the more massive fragment
    
    # Verify the values
    expected_m1c2 = 198.0
    expected_m2c2 = 99.0
    if not (math.isclose(m1c2, expected_m1c2) and math.isclose(m2c2, expected_m2c2)):
        return f"Incorrect fragment rest-mass energies. Calculated m1c2={m1c2:.2f} GeV, m2c2={m2c2:.2f} GeV. Expected 198 GeV and 99 GeV."

    # 3. Calculate total kinetic energy released (Q-value)
    Q_value_GeV = Mc2 - final_rest_mass_energy_total
    
    expected_Q_value = 3.0
    if not math.isclose(Q_value_GeV, expected_Q_value):
        return f"Incorrect Q-value. Calculated {Q_value_GeV:.2f} GeV. Expected {expected_Q_value:.2f} GeV."

    # 4. Calculate T1 using the classical (non-relativistic) approximation
    # From conservation of momentum (p1=p2), T_classical is inversely proportional to mass.
    # T1_classical / T2_classical = m2 / m1 = m2c2 / m1c2 = 1/2
    # So, T2_classical = 2 * T1_classical
    # We also know T1_classical + T2_classical = Q_value_GeV
    # T1_classical + 2 * T1_classical = Q_value_GeV
    # 3 * T1_classical = Q_value_GeV
    T1_classical_GeV = Q_value_GeV / 3.0
    
    expected_T1_classical = 1.0
    if not math.isclose(T1_classical_GeV, expected_T1_classical):
        return f"Incorrect classical T1. Calculated {T1_classical_GeV:.2f} GeV. Expected {expected_T1_classical:.2f} GeV."

    # 5. Calculate the correct (relativistic) T1
    # From conservation of momentum (p1=p2), (p1c)^2 = (p2c)^2.
    # The relativistic relation is (pc)^2 = T^2 + 2*T*(mc^2).
    # So, T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
    # We also know T2 = Q_value_GeV - T1.
    # T1^2 + 2*T1*(198) = (3 - T1)^2 + 2*(3 - T1)*(99)
    # T1^2 + 396*T1 = (9 - 6*T1 + T1^2) + (594 - 198*T1)
    # The T1^2 terms cancel.
    # 396*T1 = 603 - 204*T1
    # 600*T1 = 603
    T1_relativistic_GeV = 603.0 / 600.0
    
    expected_T1_relativistic = 1.005
    if not math.isclose(T1_relativistic_GeV, expected_T1_relativistic):
        return f"Incorrect relativistic T1. Calculated {T1_relativistic_GeV:.4f} GeV. Expected {expected_T1_relativistic:.4f} GeV."

    # 6. Find the difference and convert units
    difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    difference_MeV = difference_GeV * 1000.0
    
    # The question's options are A) 20, B) 2, C) 10, D) 5.
    # The provided answer is <<<D>>>, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0
    
    if math.isclose(difference_MeV, expected_difference_MeV):
        return "Correct"
    else:
        return (f"The calculated difference is {difference_MeV:.2f} MeV, "
                f"but the expected answer is {expected_difference_MeV} MeV.")

# Run the check
result = check_answer()
print(result)