import math

def check_fission_energy_calculation():
    """
    This function verifies the step-by-step calculation of the kinetic energy difference
    in the given nuclear fission problem.
    """
    # --- Step 1: Define given information from the problem ---
    # All energy values are in GeV for consistency in calculations.
    Mc2 = 300.0  # Initial rest-mass energy in GeV
    final_mass_fraction = 0.99
    mass_ratio_m1_to_m2 = 2.0  # m1 is the more massive fragment

    # The final answer from the LLM, in MeV. Option A is 5 MeV.
    llm_answer_value_MeV = 5.0

    # --- Step 2: Calculate the rest-mass energies of the fragments ---
    # Total rest-mass energy of fragments: (m1 + m2)c^2 = 0.99 * Mc^2
    m1c2_plus_m2c2 = final_mass_fraction * Mc2
    
    # From m1 = 2*m2, we have m1c^2 = 2*m2c^2.
    # Substituting into the sum: 2*m2c^2 + m2c^2 = m1c2_plus_m2c2
    # (mass_ratio + 1) * m2c^2 = m1c2_plus_m2c2
    m2c2 = m1c2_plus_m2c2 / (mass_ratio_m1_to_m2 + 1)
    m1c2 = mass_ratio_m1_to_m2 * m2c2

    # --- Step 3: Calculate the total kinetic energy released (Q-value) ---
    # Q = Initial Energy - Final Rest-Mass Energy
    Q = Mc2 - m1c2_plus_m2c2

    # --- Step 4: Calculate the correct (relativistic) kinetic energy T1 ---
    # From conservation of momentum, p1 = p2, so (p1c)^2 = (p2c)^2.
    # The relativistic energy-momentum relation gives (pc)^2 = T^2 + 2*T*(mc^2).
    # Therefore, T1^2 + 2*T1*(m1c^2) = T2^2 + 2*T2*(m2c^2).
    # Using T2 = Q - T1, we can solve for T1.
    # The equation simplifies to a linear equation in T1:
    # T1 * (2*m1c2 + 2*Q + 2*m2c2) = Q^2 + 2*Q*m2c2
    T1_relativistic_GeV = (Q**2 + 2*Q*m2c2) / (2*m1c2 + 2*Q + 2*m2c2)

    # --- Step 5: Calculate the classical (non-relativistic) kinetic energy T1_classical ---
    # In classical mechanics, from momentum conservation (m1*v1 = m2*v2),
    # the kinetic energies are inversely proportional to the masses: T1/T2 = m2/m1.
    # T1_cl = (m2/m1) * T2_cl
    # Q = T1_cl + T2_cl = T1_cl + (m1/m2)*T1_cl = T1_cl * (1 + m1/m2)
    T1_classical_GeV = Q / (1 + mass_ratio_m1_to_m2)

    # --- Step 6: Find the difference and compare ---
    # Convert results from GeV to MeV by multiplying by 1000
    T1_relativistic_MeV = T1_relativistic_GeV * 1000.0
    T1_classical_MeV = T1_classical_GeV * 1000.0
    
    calculated_difference_MeV = T1_relativistic_MeV - T1_classical_MeV

    # Use a small tolerance for floating-point comparison
    tolerance = 1e-6
    if abs(calculated_difference_MeV - llm_answer_value_MeV) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The final answer is incorrect.\n"
            f"The calculated difference between relativistic and classical T1 is {calculated_difference_MeV:.3f} MeV, "
            f"but the provided answer is {llm_answer_value_MeV:.3f} MeV.\n\n"
            f"Detailed calculation breakdown:\n"
            f"  - Rest-mass energy of m1 (m1c^2): {m1c2:.2f} GeV\n"
            f"  - Rest-mass energy of m2 (m2c^2): {m2c2:.2f} GeV\n"
            f"  - Total kinetic energy (Q-value): {Q:.2f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic_GeV:.4f} GeV ({T1_relativistic_MeV:.1f} MeV)\n"
            f"  - Classical T1: {T1_classical_GeV:.4f} GeV ({T1_classical_MeV:.1f} MeV)\n"
            f"  - Calculated Difference: {calculated_difference_MeV:.1f} MeV"
        )
        return reason

# Execute the check
result = check_fission_energy_calculation()
print(result)