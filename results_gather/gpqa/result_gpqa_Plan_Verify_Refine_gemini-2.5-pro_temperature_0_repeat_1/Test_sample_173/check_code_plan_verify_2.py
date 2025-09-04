import math

def check_fission_energy_calculation():
    """
    This function verifies the step-by-step calculation of the kinetic energy difference
    in a spontaneous fission problem.
    """
    # --- Problem Parameters ---
    # Initial rest-mass energy in GeV
    Mc2 = 300.0
    # Ratio of fragment rest-masses
    mass_ratio = 2.0
    # Sum of fragment rest-masses as a fraction of initial mass
    mass_sum_fraction = 0.99

    # --- LLM's Answer Details ---
    llm_final_diff_MeV = 5.0
    llm_option = 'D'

    # --- Step 1: Calculate Fragment Rest-Mass Energies ---
    # Sum of the fragments' rest-mass energies
    m1c2_plus_m2c2 = mass_sum_fraction * Mc2
    if not math.isclose(m1c2_plus_m2c2, 297.0):
        return f"Constraint check failed: Sum of fragment rest-mass energies should be 0.99 * 300 = 297 GeV, but was calculated as {m1c2_plus_m2c2} GeV."

    # From m1 = 2*m2, we have m1c^2 = 2*m2c^2.
    # Substituting into the sum: 2*m2c^2 + m2c^2 = 297 GeV => 3*m2c^2 = 297 GeV
    m2c2 = m1c2_plus_m2c2 / (mass_ratio + 1.0)
    m1c2 = mass_ratio * m2c2

    # Verify intermediate values from the LLM's response
    if not math.isclose(m1c2, 198.0):
        return f"Incorrect intermediate calculation: Rest-mass energy of the more massive fragment (m1c^2) should be 198 GeV, but was calculated as {m1c2} GeV."
    if not math.isclose(m2c2, 99.0):
        return f"Incorrect intermediate calculation: Rest-mass energy of the less massive fragment (m2c^2) should be 99 GeV, but was calculated as {m2c2} GeV."

    # --- Step 2: Calculate Total Kinetic Energy (Q-value) ---
    # The Q-value is the mass defect converted to energy.
    Q = Mc2 - m1c2_plus_m2c2
    if not math.isclose(Q, 3.0):
        return f"Incorrect intermediate calculation: Total kinetic energy (Q-value) should be 3 GeV, but was calculated as {Q} GeV."

    # --- Step 3: Calculate T1 (Relativistic) ---
    # From conservation of momentum, p1 = p2.
    # The relativistic energy-momentum relation is E^2 = (pc)^2 + (mc^2)^2.
    # This gives (pc)^2 = T^2 + 2*T*(mc^2).
    # Since (p1c)^2 = (p2c)^2:
    # T1^2 + 2*T1*m1c^2 = T2^2 + 2*T2*m2c^2
    # Using T2 = Q - T1 and solving for T1 gives:
    # T1 = (Q^2 + 2*Q*m2c^2) / (2*Mc^2)
    T1_relativistic_GeV = (Q**2 + 2 * Q * m2c2) / (2 * Mc2)
    
    if not math.isclose(T1_relativistic_GeV, 1.005):
        return f"Incorrect intermediate calculation: Relativistic T1 should be 1.005 GeV, but was calculated as {T1_relativistic_GeV} GeV."

    # --- Step 4: Calculate T1 (Classical) ---
    # In classical mechanics, T = p^2 / (2m).
    # Since p1 = p2, T1/T2 = m2/m1 = 1/2.
    # Total kinetic energy Q = T1 + T2 = T1 + 2*T1 = 3*T1.
    # Therefore, T1 = Q / 3.
    T1_classical_GeV = Q / (mass_ratio + 1.0)

    if not math.isclose(T1_classical_GeV, 1.000):
        return f"Incorrect intermediate calculation: Classical T1 should be 1.000 GeV, but was calculated as {T1_classical_GeV} GeV."

    # --- Step 5: Calculate the Difference ---
    difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    difference_MeV = difference_GeV * 1000.0

    # --- Final Verification ---
    # Use a small tolerance for floating-point comparison.
    tolerance = 1e-6
    if not math.isclose(difference_MeV, llm_final_diff_MeV, rel_tol=tolerance):
        return f"Incorrect final answer: The calculated difference is {difference_MeV:.3f} MeV, but the LLM's answer is {llm_final_diff_MeV} MeV."

    # Check if the selected option corresponds to the correct value.
    options = {'A': 2.0, 'B': 20.0, 'C': 10.0, 'D': 5.0}
    if not (llm_option in options and math.isclose(options[llm_option], llm_final_diff_MeV, rel_tol=tolerance)):
        return f"Incorrect option: The final answer is {llm_final_diff_MeV} MeV, which corresponds to option D, but the LLM selected option {llm_option}."

    return "Correct"

# Execute the check and print the result
result = check_fission_energy_calculation()
print(result)