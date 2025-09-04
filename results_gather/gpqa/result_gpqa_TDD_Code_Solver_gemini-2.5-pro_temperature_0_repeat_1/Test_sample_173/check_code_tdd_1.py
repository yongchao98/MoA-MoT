import math

def check_correctness():
    """
    This function checks the correctness of the given physics problem solution.

    The problem involves a spontaneous fission of a nucleus at rest.
    It calculates the kinetic energy of the more massive fragment using both
    relativistic and classical mechanics and finds the difference.
    """

    # --- Define problem constants from the question ---
    # Initial rest-mass energy of the nucleus in GeV
    E_M_GeV = 300.0

    # --- Step 1: Determine the rest-mass energies of the two fragments ---
    # The sum of the rest-masses of the fragments is 99% of the initial mass.
    # This also applies to their rest-mass energies.
    E_sum_fragments_GeV = 0.99 * E_M_GeV

    # One fragment (m1) is 2 times more massive than the other (m2).
    # Therefore, their rest-mass energies are in the same ratio: E_m1 = 2 * E_m2.
    # We have a system of two linear equations:
    # 1) E_m1 + E_m2 = E_sum_fragments_GeV
    # 2) E_m1 = 2 * E_m2
    #
    # Substitute (2) into (1):
    # 2 * E_m2 + E_m2 = E_sum_fragments_GeV
    # 3 * E_m2 = E_sum_fragments_GeV
    #
    # Solve for E_m2 and then E_m1.
    E_m2_GeV = E_sum_fragments_GeV / 3.0
    # E_m1 is the rest-mass energy of the more massive fragment.
    E_m1_GeV = 2.0 * E_m2_GeV

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # By conservation of energy, the initial rest-mass energy equals the sum of
    # the final rest-mass energies and the final kinetic energies.
    # E_M = (E_m1 + E_m2) + (T1 + T2)
    # The total kinetic energy Q = T1 + T2 is the "lost" mass-energy.
    Q_GeV = E_M_GeV - E_sum_fragments_GeV

    # --- Step 3: Calculate T1 using relativistic mechanics ---
    # The nucleus is initially at rest, so by conservation of momentum, the final
    # momenta of the fragments are equal in magnitude and opposite in direction.
    # |p1| = |p2| = p
    #
    # The relativistic energy-momentum relation is E^2 = (pc)^2 + (mc^2)^2.
    # Total energy E = T + mc^2 (Kinetic + Rest-mass).
    # So, (pc)^2 = (T + mc^2)^2 - (mc^2)^2 = T^2 + 2*T*mc^2.
    #
    # Apply this to both fragments:
    # (p*c)^2 = T1_rel^2 + 2 * T1_rel * E_m1
    # (p*c)^2 = T2_rel^2 + 2 * T2_rel * E_m2
    #
    # Equating these and substituting T2_rel = Q - T1_rel, we can solve for T1_rel.
    # A simplified resulting formula is:
    # T1_rel = (Q * (Q + 2*E_m2)) / (2 * E_M)
    T1_relativistic_GeV = (Q_GeV * (Q_GeV + 2.0 * E_m2_GeV)) / (2.0 * E_M_GeV)

    # --- Step 4: Calculate T1 using classical (non-relativistic) mechanics ---
    # Classical kinetic energy is T = p^2 / (2m).
    # Since momenta are equal (|p1| = |p2|), the ratio of kinetic energies is:
    # T1_cl / T2_cl = (p^2 / (2*m1)) / (p^2 / (2*m2)) = m2 / m1.
    # Since m1 = 2*m2, we have T1_cl / T2_cl = 1/2, or T2_cl = 2 * T1_cl.
    #
    # The total kinetic energy is T1_cl + T2_cl = Q.
    # Substituting T2_cl: T1_cl + 2 * T1_cl = Q => 3 * T1_cl = Q.
    # T1_cl = Q / 3.
    # A more general formula is T1_cl = Q * m2 / (m1 + m2), which in terms of energy is:
    T1_classical_GeV = Q_GeV * E_m2_GeV / (E_m1_GeV + E_m2_GeV)

    # --- Step 5: Calculate the difference and check the answer ---
    # The question asks for the difference in MeV.
    difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    calculated_difference_MeV = difference_GeV * 1000.0

    # The expected answer from option A is 5 MeV.
    expected_difference_MeV = 5.0

    # Check if the calculated value matches the expected value within a small tolerance.
    if math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-5):
        return "Correct"
    else:
        reason = (
            f"The calculated difference is incorrect.\n"
            f"The provided answer states the difference is {expected_difference_MeV} MeV.\n"
            f"My calculation results in a difference of {calculated_difference_MeV:.3f} MeV.\n\n"
            f"Intermediate values:\n"
            f"  - Initial rest-mass energy (E_M): {E_M_GeV:.2f} GeV\n"
            f"  - Sum of fragment rest-mass energies: {E_sum_fragments_GeV:.2f} GeV\n"
            f"  - Rest-mass energy of massive fragment (E_m1): {E_m1_GeV:.2f} GeV\n"
            f"  - Rest-mass energy of lighter fragment (E_m2): {E_m2_GeV:.2f} GeV\n"
            f"  - Total kinetic energy released (Q): {Q_GeV:.2f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic_GeV:.4f} GeV\n"
            f"  - Classical T1: {T1_classical_GeV:.4f} GeV\n"
            f"The calculated value of {calculated_difference_MeV:.3f} MeV does not match the expected {expected_difference_MeV} MeV."
        )
        return reason

# The final result of the check.
result = check_correctness()
print(result)