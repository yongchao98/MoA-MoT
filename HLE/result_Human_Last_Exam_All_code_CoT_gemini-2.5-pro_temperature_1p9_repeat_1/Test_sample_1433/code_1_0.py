import math

def solve_hypercube_meeting_time():
    """
    Calculates the expected time and variance for Alice and Bob to meet on a hypercube,
    and checks a related inequality.
    """
    d = 14

    # Calculate Expected Time E[X_14]
    # Let E_i be the expected time to meet from distance i.
    # We solve for Delta_i = E_i - E_{i-2}.
    # The states are i = d, d-2, ..., 2.
    delta = {}
    
    # Base case for i = d = 14
    p_d_d_minus_2 = (d * (d - 1)) / (d * d)
    delta[d] = 1 / p_d_d_minus_2

    # Solve for other delta_i recursively
    for i in range(d - 2, 0, -2):
        p_i_i_plus_2 = ((d - i) * (d - i - 1)) / (d * d)
        p_i_i_minus_2 = (i * (i - 1)) / (d * d)
        if p_i_i_minus_2 == 0:
            # This case shouldn't be reached for i > 0
            continue
        delta[i] = (1 + p_i_i_plus_2 * delta[i + 2]) / p_i_i_minus_2
    
    # E_d is the sum of all Delta_i
    expected_time = sum(delta.values())
    print(f"The integer part of the expected time EX_14 is: {int(expected_time)}")

    # Calculate Variance D^2[X_14]
    # First, we need the values of E_i for all i.
    E = {}
    E[0] = 0
    current_sum = 0
    for i in range(2, d + 1, 2):
        current_sum += delta[i]
        E[i] = current_sum

    # Let S_i = E[M_i^2]. We solve for Sigma_i = S_i - S_{i-2}.
    Sigma = {}
    
    # Base case for i = d
    p_d_d_minus_2 = (d * (d - 1)) / (d * d)
    Sigma[d] = (2 * E[d] - 1) / p_d_d_minus_2
    
    # Solve for other Sigma_i recursively
    for i in range(d - 2, 0, -2):
        p_i_i_plus_2 = ((d - i) * (d - i - 1)) / (d * d)
        p_i_i_minus_2 = (i * (i - 1)) / (d * d)
        if p_i_i_minus_2 == 0:
            continue
        Sigma[i] = (2 * E[i] - 1 + p_i_i_plus_2 * Sigma[i + 2]) / p_i_i_minus_2
    
    # S_d is the sum of all Sigma_i
    second_moment = sum(Sigma.values())
    variance = second_moment - E[d]**2
    print(f"The integer part of the variance D^2X_14 is: {int(variance)}")

    # Answer EX_15
    print("The expected time EX_15 is: infinity")
    
    # Check the inequality
    lhs = expected_time
    try:
        rhs = (d / 2) * (d**d / math.factorial(d))
    except OverflowError:
        # For large d, need to use log properties
        log_rhs = math.log(d / 2) + d * math.log(d) - math.lgamma(d + 1)
        rhs = math.exp(log_rhs)

    is_true = "yes" if lhs <= rhs else "no"
    print(f"Is it true that EX_14 <= (14/2) * (14^14 / 14!) ?")
    print(f"{int(lhs)} <= {int(rhs)} ? {is_true}")

solve_hypercube_meeting_time()