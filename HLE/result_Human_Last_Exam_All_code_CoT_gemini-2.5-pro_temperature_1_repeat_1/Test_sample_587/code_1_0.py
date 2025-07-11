import math

def combinations(n, k):
    """Calculates the binomial coefficient C(n, k)"""
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

def double_factorial(n):
    """Calculates the double factorial n!!"""
    if n < 0:
        return 0
    if n == 0 or n == 1:
        return 1
    res = 1
    for i in range(n, 1, -2):
        res *= i
    return res

def count_k_matchings_in_kn(n, k):
    """
    Calculates the number of k-matchings in a complete graph K_n.
    This is FPT for parameter k on the class of complete graphs.
    """
    print(f"Calculating the number of {k}-matchings in a complete graph K_{n}.")
    
    if 2 * k > n:
        print(f"A {k}-matching requires 2*{k}={2*k} vertices, but the graph only has {n}. Result is 0.")
        return 0
        
    # The number of ways to choose 2k vertices from n
    num_choose_vertices = combinations(n, 2 * k)
    
    # The number of ways to form a perfect matching on 2k vertices
    num_perfect_matchings = double_factorial(2 * k - 1)
    
    # Total number of k-matchings
    total_matchings = num_choose_vertices * num_perfect_matchings
    
    print("The formula is: (C(n, 2k)) * (2k-1)!!")
    print(f"Substituting the values n={n} and k={k}:")
    print(f"C({n}, {2*k}) * ({2*k - 1})!! = {num_choose_vertices} * {num_perfect_matchings} = {total_matchings}")
    
    return total_matchings

# Example demonstrating the FPT nature for the class G of cliques.
# This calculation is polynomial in n for fixed k.
n_val = 10
k_val = 3
count_k_matchings_in_kn(n_val, k_val)