import math

def check_fission_energy_calculation():
    """
    This function checks the correctness of the LLM's answer to the physics problem.
    It recalculates all values from the problem statement and compares the final result.
    """
    # --- Define constants and constraints from the problem ---
    try:
        # Initial rest-mass energy of the nucleus
        Mc2 = 300.0  # GeV

        # The sum of the rest-masses of the two fragments is 99% of the initial mass M.
        mass_sum_fraction = 0.99

        # One fragment (m1) is 2 times more massive than the other (m2).
        mass_ratio = 2.0  # m1 = 2 * m2
    except Exception as e:
        return f"Failed to define initial constants. Error: {e}"

    # --- Step 1: Calculate rest-mass energies of the fragments ---
    # (m1c^2 + m2c^2) = 0.99 * Mc^2
    m1c2_plus_m2c2 = mass_sum_fraction * Mc2
    
    # We have a system of two equations:
    # 1) m1c2 + m2c2 = m1c2_plus_m2c2
    # 2) m1c2 = mass_ratio * m2c2
    # Substitute (2) into (1):
    # (mass_ratio * m2c2) + m2c2 = m1c2_plus_m2c2
    # m2c2 * (mass_ratio + 1) = m1c2_plus_m2c2
    m2c2 = m1c2_plus_m2c2 / (mass_ratio + 1)
    m1c2 = mass_ratio * m2c2

    # --- Step 2: Calculate total kinetic energy released (Q-value) ---
    # This is the energy equivalent of the mass defect.
    T_total = Mc2 - m1c2_plus_m2c2

    # --- Step 3: Calculate the relativistic kinetic energy of the more massive fragment (T1) ---
    # The nucleus is initially at rest, so by conservation of momentum, p1 = p2.
    # The relativistic relation is (pc)^2 = E^2 - (mc^2)^2 = (T + mc^2)^2 - (mc^2)^2 = T^2 + 2*T*(mc^2).
    # So, T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # We also know from conservation of energy that T2 = T_total - T1.
    # Substituting T2 and solving for T1 yields:
    # T1 = (T_total^2 + 2*T_total*m2c2) / (2*m1c2 + 2*m2c2 + 2*T_total)
    # A simpler form is T1 = T_total * (T_total + 2*m2c2) / (2 * Mc2)
    # Let's use the LLM's derived formula for a direct check:
    T1_relativistic = (T_total**2 + 2 * T_total * m2c2) / (2 * m1c2 + 2 * m2c2 + 2 * T_total)

    # --- Step 4: Calculate the classical kinetic energy of the more massive fragment (T1) ---
    # In classical mechanics, T = p^2 / (2m), so p^2 = 2mT.
    # Since p1 = p2, we have p1^2 = p2^2, which means 2*m1*T1 = 2*m2*T2.
    # m1*T1 = m2*T2  => T2 = (m1/m2) * T1 = mass_ratio * T1
    # The total kinetic energy is T_total = T1 + T2 = T1 + (mass_ratio * T1) = T1 * (1 + mass_ratio).
    # So, T1_classical = T_total / (1 + mass_ratio).
    T1_classical = T_total / (1 + mass_ratio)

    # --- Step 5: Calculate the difference and check the answer ---
    # The question asks for the difference in MeV. 1 GeV = 1000 MeV.
    difference_MeV = (T1_relativistic - T1_classical) * 1000

    # The LLM's answer is D, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # Check if the calculated values match the LLM's intermediate steps and final answer.
    # Use a small tolerance for floating-point comparisons.
    tolerance = 1e-6

    # Check intermediate values from the LLM's response
    llm_m1c2 = 198.0
    llm_m2c2 = 99.0
    llm_T_total = 3.0
    llm_T1_relativistic = 1.0050
    llm_T1_classical = 1.0000

    if not math.isclose(m1c2, llm_m1c2, rel_tol=tolerance):
        return f"Incorrect. The rest-mass energy of the more massive fragment (m1c^2) is calculated to be {m1c2:.4f} GeV, but the LLM's value is {llm_m1c2:.4f} GeV."
    if not math.isclose(m2c2, llm_m2c2, rel_tol=tolerance):
        return f"Incorrect. The rest-mass energy of the less massive fragment (m2c^2) is calculated to be {m2c2:.4f} GeV, but the LLM's value is {llm_m2c2:.4f} GeV."
    if not math.isclose(T_total, llm_T_total, rel_tol=tolerance):
        return f"Incorrect. The total kinetic energy (T_total) is calculated to be {T_total:.4f} GeV, but the LLM's value is {llm_T_total:.4f} GeV."
    if not math.isclose(T1_relativistic, llm_T1_relativistic, rel_tol=tolerance):
        return f"Incorrect. The relativistic kinetic energy (T1) is calculated to be {T1_relativistic:.4f} GeV, but the LLM's value is {llm_T1_relativistic:.4f} GeV."
    if not math.isclose(T1_classical, llm_T1_classical, rel_tol=tolerance):
        return f"Incorrect. The classical kinetic energy (T1_classical) is calculated to be {T1_classical:.4f} GeV, but the LLM's value is {llm_T1_classical:.4f} GeV."
    
    # Final check on the difference
    if not math.isclose(difference_MeV, expected_difference_MeV, rel_tol=tolerance):
        return f"Incorrect. The final calculated difference is {difference_MeV:.4f} MeV, but the answer 'D' implies {expected_difference_MeV} MeV."

    # If all checks pass, the answer is correct.
    return "Correct"

# Run the check
result = check_fission_energy_calculation()
print(result)