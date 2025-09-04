import math

def check_fission_answer():
    """
    This function verifies the solution to the fission problem by recalculating all values
    from the initial constraints.
    """
    # --- Problem Constraints ---
    M_energy = 300.0  # Initial rest-mass energy in GeV
    mass_sum_fraction = 0.99
    m1_to_m2_ratio = 2.0

    # The answer to check is B, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # --- Step 1: Calculate rest-mass energies of the fragments ---
    # m1c^2 + m2c^2 = 0.99 * Mc^2
    m_sum_energy = mass_sum_fraction * M_energy
    
    # m1c^2 = 2 * m2c^2
    # Substitute into the sum: 2*m2c^2 + m2c^2 = m_sum_energy
    # 3 * m2c^2 = m_sum_energy
    m2_energy = m_sum_energy / (m1_to_m2_ratio + 1)
    m1_energy = m1_to_m2_ratio * m2_energy

    # --- Constraint Check: Verify fragment masses from the LLM's reasoning ---
    llm_m1_energy = 198.0
    llm_m2_energy = 99.0
    if not (math.isclose(m1_energy, llm_m1_energy) and math.isclose(m2_energy, llm_m2_energy)):
        return f"Incorrect fragment mass calculation. Expected m1={llm_m1_energy}, m2={llm_m2_energy}. Calculated m1={m1_energy:.2f}, m2={m2_energy:.2f}."

    # --- Step 2: Calculate total kinetic energy released (Q-value) ---
    # Q = T1 + T2 = Initial_Energy - Final_Rest_Energy
    Q_value = M_energy - m_sum_energy
    
    # --- Constraint Check: Verify Q-value ---
    if not math.isclose(Q_value, 3.0):
        return f"Incorrect Q-value calculation. Expected 3.0 GeV, but calculated {Q_value:.2f} GeV."

    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # Momentum conservation: p1 = p2 => sqrt(2*m1*T1) = sqrt(2*m2*T2) => m1*T1 = m2*T2
    # T1/T2 = m2/m1 = 1/2 => T2 = 2*T1
    # Energy conservation: T1 + T2 = Q => T1 + 2*T1 = Q => 3*T1 = Q
    T1_classical = Q_value / (m1_to_m2_ratio + 1)

    # --- Constraint Check: Verify classical T1 ---
    llm_T1_classical = 1.0
    if not math.isclose(T1_classical, llm_T1_classical):
        return f"Incorrect classical T1. Expected {llm_T1_classical} GeV, but calculated {T1_classical:.4f} GeV."

    # --- Step 4: Calculate T1 using the correct relativistic formula ---
    # Momentum conservation: (p1c)^2 = (p2c)^2
    # Relativistic relation: (pc)^2 = E^2 - (mc^2)^2 = (T + mc^2)^2 - (mc^2)^2 = T^2 + 2*T*mc^2
    # T1^2 + 2*T1*m1c^2 = T2^2 + 2*T2*m2c^2
    # Substitute T2 = Q - T1 and solve for T1. The general solution is:
    # T1 = (Q^2 + 2*Q*m2c^2) / (2 * (m1c^2 + m2c^2 + Q))
    # Since m1c^2 + m2c^2 + Q = Mc^2, the formula simplifies to:
    T1_relativistic = (Q_value**2 + 2 * Q_value * m2_energy) / (2 * M_energy)

    # --- Constraint Check: Verify relativistic T1 ---
    llm_T1_relativistic = 1.005
    if not math.isclose(T1_relativistic, llm_T1_relativistic):
        return f"Incorrect relativistic T1. Expected {llm_T1_relativistic} GeV, but calculated {T1_relativistic:.4f} GeV."

    # --- Step 5: Calculate the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * 1000.0

    # --- Final Verification ---
    # Check if the calculated difference matches the expected answer (5 MeV)
    if math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-4):
        return "Correct"
    else:
        return (f"Incorrect. The final calculated difference is {difference_MeV:.2f} MeV, "
                f"which does not match the expected answer of {expected_difference_MeV} MeV.")

# Run the check
result = check_fission_answer()
print(result)