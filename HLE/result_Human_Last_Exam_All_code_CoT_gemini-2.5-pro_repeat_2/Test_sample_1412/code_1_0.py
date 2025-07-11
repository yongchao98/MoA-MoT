import math

def get_prime_factorization(n):
    """
    Returns a dictionary of prime factors of n and their exponents.
    Example: get_prime_factorization(12) returns {2: 2, 3: 1}
    """
    factors = {}
    d = 2
    temp_n = n
    while d * d <= temp_n:
        while temp_n % d == 0:
            factors[d] = factors.get(d, 0) + 1
            temp_n //= d
        d += 1
    if temp_n > 1:
        factors[temp_n] = factors.get(temp_n, 0) + 1
    return factors

def num_divisors(factors):
    """
    Calculates the number of divisors from a prime factorization dictionary.
    d(p1^a1 * p2^a2) = (a1+1)*(a2+1)
    """
    if not factors:
        return 1
    res = 1
    for p in factors:
        res *= (factors[p] + 1)
    return res

def num_odd_divisors(factors):
    """
    Calculates the number of odd divisors from a prime factorization dictionary.
    This is equivalent to the number of divisors of the odd part of the number.
    """
    odd_factors = {p: exp for p, exp in factors.items() if p != 2}
    return num_divisors(odd_factors)

def solve():
    """
    Calculates the total number of non-isomorphic graphs based on the classification.
    """
    num_vertices = 2000
    # The parameter m for classification is half the number of vertices.
    m = num_vertices // 2

    # --- Type I graphs (k=2) ---
    # The number of non-isomorphic graphs is 1 if m is even, 2 if m is odd.
    num_type_I = 1 if m % 2 == 0 else 2

    # --- Type II graphs (k=0) ---
    # The number of non-isomorphic connected graphs of this type is always 1.
    num_type_II = 1

    # --- Type III graphs (k=1) ---
    # The number is given by the formula (d(m) + d_odd(m)) / 2.
    m_factors = get_prime_factorization(m)
    d_m = num_divisors(m_factors)
    d_odd_m = num_odd_divisors(m_factors)
    num_type_III = (d_m + d_odd_m) // 2

    # --- Total number of graphs ---
    total = num_type_I + num_type_II + num_type_III

    print(f"The number of vertices in the graph G is {num_vertices}.")
    print(f"The problem is solved by classifying graphs into three types based on a parameter m = {num_vertices} / 2 = {m}.")
    print("-" * 20)
    print(f"Number of Type I graphs = {num_type_I}")
    print(f"Number of Type II graphs = {num_type_II}")
    print(f"Number of Type III graphs = (d({m}) + d_odd({m})) / 2 = ({d_m} + {d_odd_m}) / 2 = {num_type_III}")
    print("-" * 20)
    print(f"Total number of non-isomorphic graphs = {num_type_I} + {num_type_II} + {num_type_III} = {total}")

solve()
<<<12>>>