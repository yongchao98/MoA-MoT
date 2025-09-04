import math

def check_physics_problem():
    """
    This function verifies the solution to the fission problem.
    It calculates the kinetic energy of the more massive fragment using both
    relativistic and classical mechanics, finds the difference, and compares
    it to the provided answer.
    """
    # --- Given values from the question ---
    # Initial rest-mass energy of the nucleus in GeV
    Mc2 = 300.0

    # --- Step 1: Calculate the rest-mass energies of the fragments ---
    # The sum of the fragments' rest-mass energies is 99% of the initial energy.
    # m1c^2 + m2c^2 = 0.99 * Mc^2
    m1c2_plus_m2c2 = 0.99 * Mc2

    # One fragment (m1) is twice as massive as the other (m2).
    # m1 = 2 * m2  =>  m1c^2 = 2 * m2c^2
    # Substitute this into the sum: 2*m2c^2 + m2c^2 = m1c2_plus_m2c2
    # 3 * m2c^2 = m1c2_plus_m2c2
    m2c2 = m1c2_plus_m2c2 / 3.0
    m1c2 = 2.0 * m2c2

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # The Q-value is the mass defect converted to energy.
    Q = Mc2 - m1c2_plus_m2c2

    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # From conservation of momentum (p1=p2) and T = p^2/(2m), we get T1/T2 = m2/m1 = 1/2.
    # So, T2_class = 2 * T1_class.
    # Since T1_class + T2_class = Q, we have T1_class + 2*T1_class = Q.
    # 3 * T1_class = Q
    T1_class = Q / 3.0

    # --- Step 4: Calculate the correct (relativistic) T1 ---
    # From conservation of momentum (p1=p2) and the relativistic relation (pc)^2 = T^2 + 2*T*mc^2,
    # we get: T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2.
    # With T2 = Q - T1, this simplifies to a linear equation for T1.
    # A more direct way is to use conservation of total energy: E1 + E2 = Mc2
    # and the momentum relation: E1^2 - (m1c2)^2 = E2^2 - (m2c2)^2
    # (E1-E2)(E1+E2) = (m1c2)^2 - (m2c2)^2
    # E1 - E2 = ((m1c2)**2 - (m2c2)**2) / Mc2
    E1_plus_E2 = Mc2
    E1_minus_E2 = (m1c2**2 - m2c2**2) / Mc2
    # Solving the system of equations for E1:
    E1 = (E1_plus_E2 + E1_minus_E2) / 2.0
    # The kinetic energy T1 is the total energy E1 minus the rest-mass energy m1c2.
    T1_rel = E1 - m1c2

    # --- Step 5: Find the difference and convert to MeV ---
    difference_GeV = T1_rel - T1_class
    difference_MeV = difference_GeV * 1000.0

    # --- Step 6: Check against the provided answer ---
    # The provided answer is <<<B>>>, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0

    # Use a small tolerance for floating-point comparisons.
    if not math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return (f"Incorrect. The calculated difference is {difference_MeV:.3f} MeV, "
                f"but the answer corresponding to <<<B>>> is {expected_difference_MeV} MeV.")

    return "Correct"

# Run the check
result = check_physics_problem()
print(result)