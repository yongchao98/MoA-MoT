import math
from fractions import Fraction

def solve_hypercube_meeting_time(d):
    """
    Calculates the expected time and variance for Alice and Bob to meet on a hypercube.
    """
    if d % 2 != 0:
        return "inf", "inf"

    # Step 1: Calculate expected times E_k
    # E_k = E[time to meet | distance is k]
    # We solve for Delta_k = E_k - E_{k-2}, where E_0 = 0.
    # E_k = sum_{j=1 to k/2} Delta_{2j}
    deltas = {}
    
    # Base case from recurrence: k=d
    deltas[d] = Fraction(d, d - 1)
    
    # Recursive calculation for Delta_k from k=d-2 down to 2
    for k in range(d - 2, 0, -2):
        numerator = d**2 + (d - k) * (d - k - 1) * deltas[k + 2]
        denominator = k * (k - 1)
        deltas[k] = Fraction(numerator, denominator)
    
    # Calculate all E_k from the deltas
    E = {}
    E[0] = Fraction(0)
    current_sum = Fraction(0)
    for k in range(2, d + 2, 2):
        current_sum += deltas[k]
        E[k] = current_sum
        
    expected_time = E[d]

    # Step 2: Calculate second moments F_k = E[T_k^2]
    # This involves solving a tridiagonal system of linear equations for F_2, F_4, ..., F_d.
    # (1-p_kk)F_k - p_k,k-2 F_{k-2} - p_k,k+2 F_{k+2} = 1 + 2(p_k,k-2 E_{k-2} + p_k,k+2 E_{k+2})
    
    N = d // 2
    M = [[Fraction(0) for _ in range(N)] for _ in range(N)]
    b = [Fraction(0) for _ in range(N)]
    d_sq = d**2

    for i in range(N):
        k = 2 * (i + 1)
        
        p_km2_num = k * (k - 1)
        p_km2 = Fraction(p_km2_num, d_sq)
        
        p_kp2_num = (d - k) * (d - k - 1)
        # Ensure numerator is non-negative (for k=d, d-k-1 is -1)
        if p_kp2_num < 0: p_kp2_num = 0
        p_kp2 = Fraction(p_kp2_num, d_sq)
        
        p_kk = Fraction(d + 2 * k * (d - k), d_sq)

        e_km2 = E.get(k - 2, Fraction(0))
        e_kp2 = E.get(k + 2, Fraction(0))
        rhs = 1 + 2 * (p_km2 * e_km2 + p_kp2 * e_kp2)
        b[i] = rhs
        
        M[i][i] = 1 - p_kk
        if i > 0:
            M[i][i-1] = -p_km2
        if i < N - 1:
            M[i][i+1] = -p_kp2

    # Step 3: Solve the tridiagonal system M*F = b using Thomas algorithm
    c_prime = [Fraction(0)] * N
    d_prime = [Fraction(0)] * N

    c_prime[0] = M[0][1] / M[0][0]
    d_prime[0] = b[0] / M[0][0]

    for i in range(1, N):
        m = M[i][i] - M[i][i-1] * c_prime[i-1]
        if i < N - 1:
            c_prime[i] = M[i][i+1] / m
        d_prime[i] = (b[i] - M[i][i-1] * d_prime[i-1]) / m

    F_sol = [Fraction(0)] * N
    F_sol[N-1] = d_prime[N-1]
    for i in range(N - 2, -1, -1):
        F_sol[i] = d_prime[i] - c_prime[i] * F_sol[i+1]
        
    F_d = F_sol[N-1]
    variance = F_d - expected_time**2
    
    return expected_time, variance

# --- Main execution ---
# Question 1: EX_14 and D^2X_14
d14 = 14
ex14, var14 = solve_hypercube_meeting_time(d14)
print(f"The integer part of the expected time EX_14 is: {int(ex14)}")
print(f"The integer part of the variance D^2X_14 is: {int(var14)}")

# Question 2: EX_15
d15 = 15
ex15, _ = solve_hypercube_meeting_time(d15)
print(f"The expected time EX_15 is: {ex15}")

# Question 3: Inequality check for d=14
d = d14
ex_d = float(ex14)
rhs = (d / 2) * (d**d / math.factorial(d))
is_true = ex_d <= rhs
print(f"Is it true that EX_d <= (d/2) * d^d/d! for d={d}? {'yes' if is_true else 'no'}")

<<<107, 7255, inf, yes>>>