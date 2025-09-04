import math

def check_correctness():
    """
    This function checks the correctness of the provided answer by recalculating the physics problem from scratch.
    
    The problem involves a spontaneous fission of a nucleus and asks for the difference
    between the relativistic and classical kinetic energy of the more massive fragment.
    """
    
    # --- Define constants and initial conditions from the problem ---
    E_M = 300.0  # Initial rest-mass energy in GeV
    mass_sum_fraction = 0.99 # The sum of fragment rest-masses is 99% of the initial mass
    mass_ratio = 2.0 # One fragment is 2 times more massive than the other
    GeV_to_MeV = 1000.0

    # --- Step 1: Calculate the rest-mass energies of the fragments ---
    # We have a system of two equations for the fragment rest-mass energies (E_m1, E_m2):
    # 1) E_m1 + E_m2 = 0.99 * E_M
    # 2) E_m1 = 2 * E_m2
    
    # Substitute (2) into (1):
    # 2 * E_m2 + E_m2 = 0.99 * E_M
    # 3 * E_m2 = 0.99 * E_M
    E_m2 = (mass_sum_fraction * E_M) / (mass_ratio + 1.0)
    E_m1 = mass_ratio * E_m2
    
    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    # The Q-value is the initial rest energy minus the final rest energy.
    Q = E_M - (E_m1 + E_m2)
    
    # --- Step 3: Calculate T1 using the classical (non-relativistic) approximation ---
    # In classical mechanics, for two fragments from a decay at rest, kinetic energy
    # is inversely proportional to mass: T1_cl / T2_cl = m2 / m1 = 1/2.
    # So, T2_cl = 2 * T1_cl.
    # Since T1_cl + T2_cl = Q, we have T1_cl + 2*T1_cl = Q, which means 3*T1_cl = Q.
    T1_classical = Q / (mass_ratio + 1.0)
    
    # --- Step 4: Calculate the correct (relativistic) T1 ---
    # From conservation of momentum (p1=p2) and the relativistic energy-momentum relation:
    # (pc)^2 = T^2 + 2*T*(mc^2)
    # We set (p1c)^2 = (p2c)^2:
    # T1^2 + 2*T1*E_m1 = T2^2 + 2*T2*E_m2
    # where T2 = Q - T1.
    # After expanding and simplifying, the T1^2 terms cancel, leaving a linear equation for T1.
    # The solution is T1 = (Q^2 + 2*Q*E_m2) / (2 * (E_m1 + E_m2 + Q))
    # A simpler form is T1 = (Q * (Q + 2*E_m2)) / (2 * E_M)
    T1_relativistic = (Q * (Q + 2 * E_m2)) / (2 * E_M)
    
    # --- Step 5: Find the difference and convert to MeV ---
    difference_GeV = T1_relativistic - T1_classical
    calculated_difference_MeV = difference_GeV * GeV_to_MeV
    
    # --- Step 6: Check against the provided answer ---
    # The provided answer is <<<A>>>, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0
    
    if not math.isclose(calculated_difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        # If the calculation does not match the answer, return the reason.
        error_message = (
            f"The provided answer 'A' (5 MeV) is incorrect.\n"
            f"The calculated difference is {calculated_difference_MeV:.3f} MeV.\n\n"
            f"Calculation details:\n"
            f"  - Rest-mass energy of massive fragment (E_m1): {E_m1:.2f} GeV\n"
            f"  - Rest-mass energy of lighter fragment (E_m2): {E_m2:.2f} GeV\n"
            f"  - Total kinetic energy released (Q): {Q:.2f} GeV\n"
            f"  - Classical T1: {T1_classical:.4f} GeV\n"
            f"  - Relativistic T1: {T1_relativistic:.4f} GeV\n"
            f"  - Difference: ({T1_relativistic:.4f} - {T1_classical:.4f}) GeV = {difference_GeV:.4f} GeV"
        )
        return error_message
        
    return "Correct"

# Run the check
result = check_correctness()
print(result)