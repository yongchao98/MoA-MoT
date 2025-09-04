import math

def check_fission_answer():
    """
    This function directly calculates the solution to the physics problem
    to verify the correctness of the provided answer.
    """
    # The provided answer is D, which corresponds to a difference of 5 MeV.
    expected_difference_MeV = 5.0

    # --- Step 1: Define initial parameters from the question (in GeV) ---
    M_c2 = 300.0

    # --- Step 2: Calculate fragment rest energies ---
    # m1 + m2 = 0.99 * M and m1 = 2 * m2
    m_fragments_sum_c2 = 0.99 * M_c2
    m2_c2 = m_fragments_sum_c2 / 3.0
    m1_c2 = 2.0 * m2_c2

    # --- Step 3: Calculate the Q-value (total kinetic energy released) ---
    Q = M_c2 - m_fragments_sum_c2

    # --- Step 4: Calculate the classical kinetic energy of the massive fragment (T1) ---
    # Classically, T1/T2 = m2/m1 = 1/2. With T1+T2=Q, this gives T1 = Q/3.
    T1_classical = Q / 3.0

    # --- Step 5: Calculate the relativistic kinetic energy of the massive fragment (T1) ---
    # From conservation of energy and momentum, we derive the exact formula for T1.
    # (p1*c)^2 = T1^2 + 2*T1*m1_c2
    # (p2*c)^2 = T2^2 + 2*T2*m2_c2
    # Since p1=p2 and T2=Q-T1, we can solve for T1.
    numerator = Q**2 + 2 * Q * m2_c2
    denominator = 2 * (m1_c2 + m2_c2) + 2 * Q
    T1_relativistic = numerator / denominator

    # --- Step 6: Calculate the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    calculated_difference_MeV = difference_GeV * 1000.0

    # --- Step 7: Check if the calculated difference matches the expected answer ---
    if math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-5):
        return "Correct"
    else:
        reason = (
            f"The answer is incorrect. The expected difference is {expected_difference_MeV} MeV.\n"
            f"The calculated values are:\n"
            f" - Classical T1 = {T1_classical:.4f} GeV\n"
            f" - Relativistic T1 = {T1_relativistic:.4f} GeV\n"
            f" - Calculated difference = ({T1_relativistic:.4f} - {T1_classical:.4f}) GeV * 1000 = {calculated_difference_MeV:.2f} MeV.\n"
            f"The calculated difference of {calculated_difference_MeV:.2f} MeV does not match the expected {expected_difference_MeV} MeV."
        )
        return reason

# Run the check
result = check_fission_answer()
print(result)