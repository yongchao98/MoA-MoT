import numpy as np
import math

def solve_hypercube_problem():
    """
    Solves the hypercube meeting time problem for d=14 and d=15.
    """
    
    # Part 1: Expected time for d=14
    d14 = 14
    ex14 = 0
    # The expected time EX_d for even d is the sum of 2j/(2j-1) for j from 1 to d/2
    # EX_14 = 2/1 + 4/3 + 6/5 + ... + 14/13
    ex_k_vals = {}
    current_ex_sum = 0
    for j in range(1, (d14 // 2) + 1):
        k = 2 * j
        current_ex_sum += k / (k - 1)
        ex_k_vals[k] = current_ex_sum
    ex14 = ex_k_vals[d14]
    
    print(f"Integer part of EX_14: {math.floor(ex14)}")

    # Part 2: Variance for d=14
    # We solve a system of linear equations for F_k = E[T^2 | dist=k]
    N = d14 // 2
    A = np.zeros((N, N))
    b = np.zeros(N)

    # Setup the tridiagonal system A*F = b
    for i in range(N):
        k = 2 * (i + 1)
        
        # Coefficients based on transition probabilities
        d_k = k * (k - 1) + (d14 - k) * (d14 - k - 1)
        p1_k = k * (k - 1)
        p2_k = (d14 - k) * (d14 - k - 1)

        A[i, i] = d_k
        if i > 0:
            A[i, i - 1] = -p1_k
        if i < N - 1:
            A[i, i + 1] = -p2_k
        
        b[i] = d14**2 * (1 + 2 * ex_k_vals[k])

    # Solve for F_k values and compute variance
    F_vec = np.linalg.solve(A, b)
    f14 = F_vec[-1]
    var14 = f14 - ex14**2
    
    print(f"Integer part of D^2(X_14): {math.floor(var14)}")

    # Part 3: Expected time for d=15
    # For odd dimensions, the parity of the distance is invariant.
    # Initial distance is 15 (odd), meeting distance is 0 (even). They never meet.
    print("EX_15: infinity")

    # Part 4: Inequality check for d=14
    # Is EX_d <= (d/2) * d^d / d! ?
    lhs = ex14
    
    # Use logarithms to handle large numbers in the factorial and power
    log_rhs = math.log(d14 / 2) + d14 * math.log(d14) - math.lgamma(d14 + 1)
    rhs = math.exp(log_rhs)
    
    is_true = "yes" if lhs <= rhs else "no"
    print(f"Is EX_14 <= 14/2 * 14^14 / 14! ? {is_true}")

solve_hypercube_problem()