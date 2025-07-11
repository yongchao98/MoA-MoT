import math

def factorial(n):
    """Computes n!"""
    return 1 if n == 0 else math.factorial(n)

def combinations(n, k):
    """Computes nCk"""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def num_involutions(k):
    """
    Computes the number of elements x in S_k such that x^2 = 1.
    These are permutations whose cycle decomposition contains only 1-cycles and 2-cycles.
    The number follows the recurrence a_k = a_{k-1} + (k-1)*a_{k-2}.
    """
    if k == 0: return 1
    a = [0] * (k + 1)
    a[0] = 1
    for i in range(1, k + 1):
        val = a[i-1]
        if i >= 2:
            val += (i-1) * a[i-2]
        a[i] = val
    return a[k]

def num_order_div_5(k):
    """
    Computes the number of elements y in S_k such that y^5 = 1.
    For k < 10, such an element is either the identity or a 5-cycle.
    """
    if k < 5:
        return 1
    # Number of 5-cycles is C(k, 5) * (5-1)!
    num_5_cycles = combinations(k, 5) * factorial(4)
    return 1 + num_5_cycles

def solve_subgroup_problem():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    print(f"Finding the number of subgroups of index {N} in G = C_2 * C_5.\n")

    a = [num_involutions(k) for k in range(N + 1)]
    b = [num_order_div_5(k) for k in range(N + 1)]
    h = [a_k * b_k for a_k, b_k in zip(a, b)]

    print("Step 1: Calculate h_k = |Hom(G, S_k)| for k=1 to 7")
    for k in range(1, N + 1):
        print(f"h_{k} = a_{k} * b_{k} = {a[k]} * {b[k]} = {h[k]}")
    print("-" * 20)

    t = [0] * (N + 1)
    h[0] = 1  # Base case for the recurrence

    print("Step 2: Calculate t_k = number of transitive homomorphisms for k=1 to 7")
    for n in range(1, N + 1):
        sum_val = 0
        for k in range(1, n):
            term = combinations(n - 1, k - 1) * t[k] * h[n-k]
            sum_val += term
        t[n] = h[n] - sum_val
        print(f"t_{n} = h_{n} - (sum) = {h[n]} - {sum_val} = {t[n]}")
    print("-" * 20)

    print(f"Step 3: Calculate the number of subgroups of index {N}")
    if t[N] < 0:
        print("Error in calculation, t_N is negative.")
        return

    num_subgroups = t[N] // factorial(N - 1)
    print(f"The number of subgroups is t_{N} / ({N}-1)!")
    print(f"Equation: {t[N]} / {factorial(N-1)} = {num_subgroups}")
    
    print("\nFinal Answer:")
    print(f"The number of subgroups of index 7 in G is {num_subgroups}.")


solve_subgroup_problem()