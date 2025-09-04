import math

def check_correctness():
    """
    This function checks the correctness of the provided answer to the physics problem.
    It calculates the kinetic energy of the more massive fragment (T1) using both
    relativistic and classical mechanics, finds the difference, and compares it
    to the value given in the selected option.
    """
    
    # --- Problem Parameters ---
    # Initial rest-mass energy of the nucleus in GeV
    E_M_GeV = 300.0
    
    # Sum of rest-masses of the two fragments is 99% of the initial mass
    mass_sum_fraction = 0.99
    
    # One fragment is 2 times more massive than the other (m1 = 2*m2)
    mass_ratio_m1_to_m2 = 2.0
    
    # The provided answer is B, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # --- Calculations ---

    # 1. Calculate the rest-mass energies of the fragments
    # Total rest-mass energy of the two fragments
    E_fragments_sum_GeV = E_M_GeV * mass_sum_fraction
    
    # Let E1 = m1*c^2 and E2 = m2*c^2.
    # We have E1 + E2 = E_fragments_sum_GeV and E1 = 2 * E2.
    # So, 3 * E2 = E_fragments_sum_GeV
    E2_GeV = E_fragments_sum_GeV / (mass_ratio_m1_to_m2 + 1)  # Rest-mass energy of the less massive fragment (m2)
    E1_GeV = mass_ratio_m1_to_m2 * E2_GeV  # Rest-mass energy of the more massive fragment (m1)

    # 2. Calculate the total kinetic energy released (Q-value)
    # This is the energy deficit from the rest masses.
    Q_GeV = E_M_GeV - E_fragments_sum_GeV

    # 3. Calculate T1 using the classical (non-relativistic) approximation
    # From classical conservation of momentum (p1 = p2): T1/T2 = m2/m1 = 1/2.
    # With T1 + T2 = Q, we get T1 + 2*T1 = Q => 3*T1 = Q.
    T1_classical_GeV = Q_GeV / 3.0

    # 4. Calculate the correct T1 using relativistic mechanics
    # From relativistic conservation of momentum (|p1| = |p2|):
    # The energy-momentum relation is E^2 = (pc)^2 + (mc^2)^2.
    # This gives (pc)^2 = (T + mc^2)^2 - (mc^2)^2 = T^2 + 2*T*(mc^2).
    # Equating (p1c)^2 and (p2c)^2:
    # T1^2 + 2*T1*E1 = T2^2 + 2*T2*E2
    # Solving for T1 with T2 = Q - T1 gives the formula:
    # T1 = (Q^2 + 2*Q*E2) / (2*E_M)
    T1_relativistic_GeV = (Q_GeV**2 + 2 * Q_GeV * E2_GeV) / (2 * E_M_GeV)

    # 5. Calculate the difference and convert to MeV
    difference_GeV = T1_relativistic_GeV - T1_classical_GeV
    calculated_difference_MeV = difference_GeV * 1000

    # 6. Compare the calculated result with the expected answer
    # Use a small tolerance for floating-point comparisons
    tolerance = 1e-3 # 0.001 MeV is a reasonable tolerance
    if abs(calculated_difference_MeV - expected_difference_MeV) < tolerance:
        return "Correct"
    else:
        reason = (
            f"The provided answer is incorrect. The calculated difference is {calculated_difference_MeV:.3f} MeV, but the answer B suggests it should be {expected_difference_MeV:.1f} MeV.\n"
            f"Calculation details:\n"
            f"  - Q-value = {Q_GeV:.1f} GeV\n"
            f"  - Rest-mass energy of lighter fragment (m2*c^2) = {E2_GeV:.1f} GeV\n"
            f"  - Classical T1 = {T1_classical_GeV:.4f} GeV\n"
            f"  - Relativistic T1 = {T1_relativistic_GeV:.4f} GeV\n"
            f"  - Difference = ({T1_relativistic_GeV:.4f} - {T1_classical_GeV:.4f}) * 1000 = {calculated_difference_MeV:.3f} MeV"
        )
        return reason

# The code block above defines the checking function.
# To get the result, we would call the function.
# result = check_correctness()
# print(result)