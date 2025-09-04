import math

def check_correctness():
    """
    This function checks the correctness of the final answer by recalculating the physics problem from scratch.
    """
    # --- Step 1: Define initial conditions and calculate fragment properties ---
    E_M = 300.0  # Initial rest-mass energy in GeV
    
    # The sum of rest-mass energies of the two fragments is 99% of the initial rest-mass energy.
    sum_final_rest_energies = 0.99 * E_M
    
    # One fragment is 2 times more massive than the other.
    # Let m1c2 and m2c2 be the rest-mass energies of the more massive and lighter fragments.
    # m1c2 = 2 * m2c2
    # m1c2 + m2c2 = sum_final_rest_energies
    # Solving this system: 3 * m2c2 = sum_final_rest_energies
    m2c2 = sum_final_rest_energies / 3.0
    m1c2 = 2.0 * m2c2

    # --- Step 2: Calculate the total kinetic energy released (Q-value) ---
    Q = E_M - sum_final_rest_energies

    # --- Step 3: Calculate T1 using Classical Approximation ---
    # From conservation of momentum (p1=p2) and T = p^2/(2m), we get T1/T2 = m2/m1 = 1/2.
    # So, T2_cl = 2 * T1_cl.
    # T1_cl + T2_cl = Q  =>  T1_cl + 2 * T1_cl = Q  =>  3 * T1_cl = Q
    T1_classical = Q / 3.0

    # --- Step 4: Calculate T1 using Relativistic Mechanics ---
    # From conservation of momentum (p1=p2) and (pc)^2 = T^2 + 2*T*(mc^2), we get:
    # T1^2 + 2*T1*(m1c2) = T2^2 + 2*T2*(m2c2)
    # We know T2 = Q - T1.
    # T1^2 + 2*T1*m1c2 = (Q - T1)^2 + 2*(Q - T1)*m2c2
    # T1^2 + 2*T1*m1c2 = Q^2 - 2*Q*T1 + T1^2 + 2*Q*m2c2 - 2*T1*m2c2
    # The T1^2 terms cancel. Rearranging to solve for T1:
    # T1 * (2*m1c2 + 2*Q + 2*m2c2) = Q^2 + 2*Q*m2c2
    # T1 * (2*(m1c2 + m2c2) + 2*Q) = Q * (Q + 2*m2c2)
    # Since (m1c2 + m2c2) + Q = E_M, the term in the parenthesis is 2*E_M.
    # T1 * (2 * E_M) = Q * (Q + 2*m2c2)
    T1_relativistic = (Q * (Q + 2 * m2c2)) / (2 * E_M)

    # --- Step 5: Calculate the difference ---
    difference_GeV = T1_relativistic - T1_classical
    
    # --- Step 6: Convert to MeV ---
    difference_MeV = difference_GeV * 1000.0
    
    # The final answer provided is <<<C>>>, which corresponds to 5 MeV.
    expected_difference_MeV = 5.0
    
    # Check if the calculated difference matches the expected value for option C.
    if math.isclose(difference_MeV, expected_difference_MeV, rel_tol=1e-6):
        return "Correct"
    else:
        return f"The calculated difference is {difference_MeV:.3f} MeV, which does not match the expected 5 MeV from option C. The provided answer is incorrect."

# The final answer from the analysis is <<<C>>>.
# The code will verify if the calculation leading to 5 MeV (Option C) is correct.
print(check_correctness())