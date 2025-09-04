import math

def check_fission_kinetic_energy():
    """
    This function verifies the calculation for the difference between relativistic and
    classical kinetic energy of a fission fragment.
    """
    # --- Step 1: Define constants from the problem statement ---
    Mc2 = 300.0  # GeV, initial rest-mass energy
    mass_sum_fraction = 0.99
    mass_ratio_m1_to_m2 = 2.0

    # --- Step 2: Calculate derived parameters based on the problem ---
    # Total rest-mass energy of the two fragments
    m1c2_plus_m2c2 = mass_sum_fraction * Mc2
    
    # Check if the sum is correct
    if not math.isclose(m1c2_plus_m2c2, 297.0):
        return f"Constraint check failed: Sum of fragment rest-mass energies should be 0.99 * 300 = 297 GeV, but is {m1c2_plus_m2c2}."

    # Individual rest-mass energies of fragments
    # We have a system of two equations:
    # 1) m1c2 + m2c2 = 297
    # 2) m1c2 = 2 * m2c2
    # Substituting (2) into (1): 2*m2c2 + m2c2 = 297 => 3*m2c2 = 297
    m2c2 = m1c2_plus_m2c2 / (mass_ratio_m1_to_m2 + 1)
    m1c2 = mass_ratio_m1_to_m2 * m2c2
    
    # Check if individual masses are correct
    if not math.isclose(m1c2, 198.0) or not math.isclose(m2c2, 99.0):
        return f"Constraint check failed: Fragment rest-mass energies should be 198 GeV and 99 GeV, but were calculated as {m1c2} and {m2c2}."

    # Total kinetic energy released
    T_total = Mc2 - m1c2_plus_m2c2
    
    if not math.isclose(T_total, 3.0):
        return f"Constraint check failed: Total kinetic energy should be 300 - 297 = 3 GeV, but was calculated as {T_total}."

    # --- Step 3: Calculate the correct (relativistic) kinetic energy for the more massive fragment (T1) ---
    # The formula derived from conservation of energy and momentum is:
    # T1 = (T_total^2 + 2 * T_total * m2c2) / (2 * Mc2)
    T1_relativistic = (T_total**2 + 2 * T_total * m2c2) / (2 * Mc2)

    # --- Step 4: Calculate the classical kinetic energy for the same fragment (T1) ---
    # The formula derived from classical conservation of momentum is:
    # T1 = T_total / (1 + m1/m2)
    T1_classical = T_total / (1 + m1c2 / m2c2)

    # --- Step 5: Calculate the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * 1000

    # --- Step 6: Verify the final answer ---
    # The LLM's code calculates a difference of 5.0 MeV, which corresponds to option D.
    expected_difference_MeV = 5.0
    
    # Check if the intermediate calculations from the LLM's code are correct
    llm_T1_relativistic = 1.005
    llm_T1_classical = 1.0
    
    if not math.isclose(T1_relativistic, llm_T1_relativistic, rel_tol=1e-6):
        return f"Calculation Error: The re-calculated relativistic T1 is {T1_relativistic:.6f} GeV, which does not match the LLM's value of {llm_T1_relativistic} GeV."
        
    if not math.isclose(T1_classical, llm_T1_classical, rel_tol=1e-6):
        return f"Calculation Error: The re-calculated classical T1 is {T1_classical:.6f} GeV, which does not match the LLM's value of {llm_T1_classical} GeV."

    # Check if the final answer is correct
    if math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-3):
        return "Correct"
    else:
        return f"Incorrect: The calculated difference is {difference_MeV:.2f} MeV, which does not match the expected answer of {expected_difference_MeV} MeV."

# Execute the check
result = check_fission_kinetic_energy()
print(result)