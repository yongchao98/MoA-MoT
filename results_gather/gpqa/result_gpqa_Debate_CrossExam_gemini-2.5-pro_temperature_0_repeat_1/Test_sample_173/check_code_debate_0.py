import math

def check_fission_energy_calculation():
    """
    This function verifies the step-by-step calculation provided in the LLM's answer.
    It recalculates all values based on the problem's constraints and compares them
    to the results in the provided solution.
    """
    # --- Define problem constants and parameters ---
    # All energy values are in GeV unless specified otherwise.
    Mc2 = 300.0  # Initial rest-mass energy
    final_mass_fraction = 0.99
    # m1 is the more massive fragment, m2 is the less massive one.
    # m1 = 2 * m2, so m1c^2 = 2 * m2c^2
    mass_ratio_m1_over_m2 = 2.0
    
    # --- Values from the LLM's answer to be checked ---
    expected_m1c2 = 198.0
    expected_m2c2 = 99.0
    expected_Q = 3.0
    expected_T1_rel = 1.005
    expected_T1_class = 1.0
    expected_difference_MeV = 5.0

    # --- Step 1: Determine Rest-Mass Energies of Fragments ---
    total_fragment_mass_energy = Mc2 * final_mass_fraction
    # We have a system of two equations:
    # 1) m1c2 + m2c2 = total_fragment_mass_energy
    # 2) m1c2 = mass_ratio_m1_over_m2 * m2c2
    # Substitute (2) into (1):
    # (mass_ratio_m1_over_m2 * m2c2) + m2c2 = total_fragment_mass_energy
    # m2c2 * (mass_ratio_m1_over_m2 + 1) = total_fragment_mass_energy
    m2c2_calc = total_fragment_mass_energy / (mass_ratio_m1_over_m2 + 1)
    m1c2_calc = mass_ratio_m1_over_m2 * m2c2_calc

    if not (math.isclose(m1c2_calc, expected_m1c2) and math.isclose(m2c2_calc, expected_m2c2)):
        return f"Error in Step 1 (Fragment Masses). Calculated m1c^2={m1c2_calc:.3f} GeV, m2c^2={m2c2_calc:.3f} GeV. Expected m1c^2={expected_m1c2} GeV, m2c^2={expected_m2c2} GeV."

    # --- Step 2: Calculate Total Kinetic Energy (Q-value) ---
    Q_calc = Mc2 - total_fragment_mass_energy
    if not math.isclose(Q_calc, expected_Q):
        return f"Error in Step 2 (Q-value). Calculated Q={Q_calc:.3f} GeV. Expected Q={expected_Q} GeV."

    # --- Step 3: Calculate Relativistic Kinetic Energy (T1_rel) ---
    # From conservation of momentum, p1 = p2, so (p1*c)^2 = (p2*c)^2.
    # Using the relation (p*c)^2 = T * (T + 2*m*c^2):
    # T1 * (T1 + 2*m1c^2) = T2 * (T2 + 2*m2c^2)
    # We also know T1 + T2 = Q, so T2 = Q - T1.
    # T1 * (T1 + 2*m1c^2) = (Q - T1) * ((Q - T1) + 2*m2c^2)
    # T1^2 + 2*m1c^2*T1 = Q^2 - Q*T1 + 2*Q*m2c^2 - Q*T1 + T1^2 - 2*m2c^2*T1
    # 2*m1c^2*T1 = Q^2 + 2*Q*m2c^2 - 2*Q*T1 - 2*m2c^2*T1
    # (2*m1c^2 + 2*Q + 2*m2c^2) * T1 = Q^2 + 2*Q*m2c^2
    # T1 = (Q^2 + 2*Q*m2c^2) / (2*m1c^2 + 2*m2c^2 + 2*Q)
    T1_rel_calc = (Q_calc**2 + 2 * Q_calc * m2c2_calc) / (2 * m1c2_calc + 2 * m2c2_calc + 2 * Q_calc)
    
    if not math.isclose(T1_rel_calc, expected_T1_rel):
        return f"Error in Step 3 (Relativistic T1). Calculated T1_rel={T1_rel_calc:.4f} GeV. Expected T1_rel={expected_T1_rel} GeV."

    # --- Step 4: Calculate Classical Kinetic Energy (T1_class) ---
    # In classical mechanics, T = p^2 / (2m).
    # Since p1 = p2, we have 2*m1*T1 = 2*m2*T2, which simplifies to m1*T1 = m2*T2.
    # Therefore, T1/T2 = m2/m1.
    # We also know T1 + T2 = Q. So, T2 = Q - T1.
    # T1 / (Q - T1) = m2 / m1
    # T1*m1 = m2*(Q - T1) = m2*Q - m2*T1
    # T1*(m1 + m2) = m2*Q
    # T1 = Q * (m2 / (m1 + m2))
    # This is equivalent to T1 = Q / (1 + m1/m2)
    T1_class_calc = Q_calc / (1 + mass_ratio_m1_over_m2)

    if not math.isclose(T1_class_calc, expected_T1_class):
        return f"Error in Step 4 (Classical T1). Calculated T1_class={T1_class_calc:.4f} GeV. Expected T1_class={expected_T1_class} GeV."

    # --- Step 5: Find the Difference ---
    difference_GeV = T1_rel_calc - T1_class_calc
    difference_MeV = difference_GeV * 1000.0

    if not math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return f"Error in Step 5 (Difference). Calculated difference={difference_MeV:.3f} MeV. Expected difference={expected_difference_MeV} MeV."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_fission_energy_calculation()
print(result)