import math

def check_fission_energy_difference():
    """
    This function verifies the calculation for the difference between relativistic and 
    classical kinetic energy of a fission fragment.

    The problem states:
    - Initial nucleus rest-mass energy (Mc^2) = 300 GeV.
    - It splits into two fragments, m1 and m2.
    - m1 = 2 * m2 (m1 is the more massive fragment).
    - m1 + m2 = 0.99 * M (sum of rest masses).
    - T1 is the kinetic energy of the more massive fragment.
    - The goal is to find the difference between the correct (relativistic) T1 
      and the classical T1.
    """

    # --- Step 1: Define initial conditions and calculate fragment rest-mass energies ---
    Mc2 = 300.0  # Initial rest-mass energy in GeV

    # From m1 + m2 = 0.99M and m1 = 2m2, we can find the rest-mass energies.
    # m1c^2 + m2c^2 = 0.99 * Mc^2
    # m1c^2 = 2 * m2c^2
    # Substituting: 2*m2c^2 + m2c^2 = 0.99 * Mc^2 => 3*m2c^2 = 0.99 * Mc^2
    m2c2 = (0.99 * Mc2) / 3.0  # Rest-mass energy of the lighter fragment
    m1c2 = 2.0 * m2c2          # Rest-mass energy of the more massive fragment

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # Q is the mass defect converted to energy.
    Q = Mc2 - (m1c2 + m2c2)

    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # In classical mechanics, from conservation of momentum (p1=p2), we have:
    # T1/T2 = m2/m1 = 1/2 => T2 = 2*T1
    # Since T1 + T2 = Q, we have T1 + 2*T1 = Q => 3*T1 = Q
    T1_classical = Q / 3.0

    # --- Step 4: Calculate the correct (relativistic) T1 ---
    # From relativistic momentum conservation (p1=p2), we use (pc)^2 = T^2 + 2*T*(mc^2).
    # T1^2 + 2*T1*m1c2 = T2^2 + 2*T2*m2c2
    # Substitute T2 = Q - T1 and solve for T1. The simplified result is:
    # T1 = (Q * (Q + 2*m2c2)) / (2 * Mc2)
    T1_relativistic = (Q * (Q + 2.0 * m2c2)) / (2.0 * Mc2)

    # --- Step 5: Find the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    difference_MeV = difference_GeV * 1000.0

    # --- Step 6: Check against the provided answer ---
    # The provided answer is <<<C>>>, which corresponds to 5 MeV from the options.
    expected_answer_value = 5.0

    # Check if the calculated value matches the expected answer from the LLM.
    # Use a small tolerance for floating-point comparison.
    if not math.isclose(difference_MeV, expected_answer_value, rel_tol=1e-6):
        return (f"Incorrect. The calculated difference is {difference_MeV:.4f} MeV, "
                f"which does not match the expected answer of {expected_answer_value} MeV "
                f"corresponding to option C.")

    # All calculations match the provided correct answer.
    return "Correct"

# Run the check
result = check_fission_energy_difference()
print(result)