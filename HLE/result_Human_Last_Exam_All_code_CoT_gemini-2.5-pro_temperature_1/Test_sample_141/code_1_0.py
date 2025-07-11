import math

def solve():
    """
    This function calculates the value of a_21 based on the described combinatorial problem.
    It uses recurrence relations for the Fubini numbers (f_n) and the main sequence (a_n).
    """
    N = 21

    # Step 1: Compute Fubini numbers f_n up to f_{N+1} = f_22,
    # as the recurrence for a_21 needs f_1 through f_21.
    # The recurrence is f_n = sum_{k=0}^{n-1} C(n,k) * f_k for n >= 1, with f_0 = 1.
    f = [0] * (N + 2)
    f[0] = 1
    for n in range(1, N + 2):
        sum_val = 0
        for k in range(n):
            sum_val += math.comb(n, k) * f[k]
        f[n] = sum_val

    # Step 2: Compute a_n up to a_N = a_21.
    # The recurrence is a_{n+1} = sum_{k=0}^{n} C(n,k) * f_{k+1} * a_{n-k}, with a_0 = 1.
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(N):  # n runs from 0 to 20
        sum_val = 0
        for k in range(n + 1):
            term = math.comb(n, k) * f[k + 1] * a[n - k]
            sum_val += term
        a[n + 1] = sum_val
    
    # Step 3: Output the final result.
    result = a[N]
    # The final equation is a_21 = <result>. The numbers in it are 21 and the result.
    print(f"a_{N} = {result}")

solve()