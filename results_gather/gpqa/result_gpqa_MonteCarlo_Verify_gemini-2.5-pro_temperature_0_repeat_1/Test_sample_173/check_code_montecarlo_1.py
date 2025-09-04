import math

def check_fission_problem_correctness():
    """
    This function checks the correctness of the given answer to the physics problem.
    The problem involves a spontaneous fission of a nucleus and asks for the difference
    between the relativistic and classical kinetic energy of the more massive fragment.
    The provided answer 'D' corresponds to a difference of 5 MeV. This code verifies that calculation.
    """

    # --- Problem Parameters from the question ---
    # Using GeV as the base unit for energy for simplicity.
    M_c2 = 300.0  # Initial rest-mass energy of the nucleus in GeV.

    # --- Step 1: Calculate the rest energies of the fragments ---
    # Let the fragments have rest energies m1_c2 and m2_c2.
    # Constraint 1: The sum of rest-masses is 99% of the initial mass.
    sum_fragments_rest_energy = 0.99 * M_c2
    
    # Constraint 2: One fragment is 2 times more massive than the other.
    # Let m1 be the more massive fragment, so m1 = 2 * m2.
    # This implies m1_c2 = 2 * m2_c2.

    # We solve the system of two equations:
    # 1) m1_c2 + m2_c2 = sum_fragments_rest_energy
    # 2) m1_c2 = 2 * m2_c2
    # Substitute (2) into (1):
    # 2 * m2_c2 + m2_c2 = sum_fragments_rest_energy
    # 3 * m2_c2 = 0.99 * M_c2
    m2_c2 = (0.99 * M_c2) / 3.0  # Rest energy of the lighter fragment.
    m1_c2 = 2.0 * m2_c2          # Rest energy of the more massive fragment.

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # By conservation of energy, the initial rest energy is converted into the
    # final rest energy plus the kinetic energy of the fragments.
    # Q = T1 + T2 = M_c2 - (m1_c2 + m2_c2)
    Q = M_c2 - sum_fragments_rest_energy
    
    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # In classical mechanics, kinetic energy T = p^2 / (2m).
    # The nucleus is initially at rest, so by conservation of momentum, the fragments
    # have equal and opposite momenta: |p1| = |p2|.
    # Therefore, p1^2 = p2^2, which means 2*m1*T1_cl = 2*m2*T2_cl.
    # m1 * T1_cl = m2 * T2_cl
    # Since m1 = 2*m2, we have T2_cl = 2 * T1_cl.
    # We also know T1_cl + T2_cl = Q.
    # Substituting T2_cl: T1_cl + 2*T1_cl = Q => 3*T1_cl = Q.
    T1_classical = Q / 3.0

    # --- Step 4: Calculate the correct T1 using relativistic mechanics ---
    # The relativistic relation is (pc)^2 = E^2 - (mc^2)^2 = (T + mc^2)^2 - (mc^2)^2.
    # This simplifies to (pc)^2 = T^2 + 2*T*(mc^2).
    # From conservation of momentum, |p1| = |p2|, so (p1c)^2 = (p2c)^2.
    # T1_rel^2 + 2*T1_rel*(m1_c2) = T2_rel^2 + 2*T2_rel*(m2_c2)
    # We also know T2_rel = Q - T1_rel.
    # Substituting T2_rel and solving for T1_rel gives:
    # T1_rel = (Q * (Q + 2 * m2_c2)) / (2 * M_c2)
    T1_relativistic = (Q * (Q + 2 * m2_c2)) / (2 * M_c2)

    # --- Step 5: Calculate the difference and check the answer ---
    # The question asks for the difference in MeV (1 GeV = 1000 MeV).
    difference_GeV = T1_relativistic - T1_classical
    calculated_difference_MeV = difference_GeV * 1000.0

    # The provided answer is 'D', which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # Check if the calculated value matches the expected value from the answer.
    if not math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        reason = (
            f"The answer 'D' (5 MeV) is incorrect.\n"
            f"The calculated difference between the relativistic and classical kinetic energy is {calculated_difference_MeV:.4f} MeV.\n"
            f"Calculation details:\n"
            f"  - Q-value (total kinetic energy): {Q:.4f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic*1000:.2f} MeV\n"
            f"  - Classical T1: {T1_classical*1000:.2f} MeV\n"
            f"  - Calculated Difference: {calculated_difference_MeV:.2f} MeV"
        )
        return reason

    return "Correct"

# Run the check and print the result.
# This will return "Correct" if the calculation confirms that 5 MeV is the right answer,
# and an error message otherwise.
print(check_fission_problem_correctness())