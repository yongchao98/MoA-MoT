import math

# Using a cache for the recursive function f(k, n) and its components
memo_cache = {}

def nCr_safe(n, r):
    """
    Safe combination function that returns 0 if r < 0 or r > n.
    Standard math.comb handles large integers perfectly.
    """
    if r < 0 or r > n:
        return 0
    return math.comb(n, r)

def A(k, n):
    """
    Calculates A_k, the number of k-subsets of the vector space (F_2)^n
    that sum to the zero vector.
    """
    if ('A', k, n) in memo_cache:
        return memo_cache[('A', k, n)]
    
    N = 2**n
    if k < 0:
        return 0
    
    # The formulas for A_k give an integer result, so integer division is safe.
    if k % 2 == 1:  # k is odd
        numerator = nCr_safe(N, k)
    else:  # k is even
        term1 = nCr_safe(N, k)
        # The combination term is non-zero only if k/2 is an integer
        term2 = (N - 1) * nCr_safe(N // 2, k // 2)
        
        # This handles the (-1)^(k/2) factor
        if (k // 2) % 2 == 1:
            numerator = term1 - term2
        else:
            numerator = term1 + term2
            
    result = numerator // N
    memo_cache[('A', k, n)] = result
    return result

def f(k, n):
    """
    Calculates f_k, the number of k-sets T satisfying the problem conditions.
    This is calculated recursively with memoization using the relation:
    f(k) = A(k) - f(k-1)
    """
    if ('f', k, n) in memo_cache:
        return memo_cache[('f', k, n)]
    
    if k == 0:
        return 1
    if k < 0:
        return 0
        
    # The recursive step
    result = A(k, n) - f(k - 1, n)
    
    # Store the result in the cache
    memo_cache[('f', k, n)] = result
    return result

def solve_and_print(n, m):
    """
    Main function to solve the problem for given n and m, and print the steps.
    """
    # Clear cache for a fresh run with new n, m values
    global memo_cache
    memo_cache.clear()

    # The final answer is f(m, n). Calling this populates the cache
    # with all necessary intermediate values of f(k,n) and A(k,n).
    final_answer = f(m, n)

    print(f"Solving for n = {n}, m = {m}\n")
    print("Let f(k) be the number of sets T of size k satisfying the conditions.")
    print("The recurrence relation is f(k) = A(k) - f(k-1).")
    print("A(k) is the number of k-subsets of the full powerset P(S) whose elements' symmetric difference is the empty set.")
    
    # Base case
    f_prev = 1
    print(f"\nBase case: f(0) = {f_prev}")
    
    # Iteratively print calculations from k=1 to m
    for k in range(1, m + 1):
        # Retrieve values from cache
        a_k = memo_cache[('A', k, n)]
        current_f = memo_cache[('f', k, n)]
        
        print(f"\nCalculating for k = {k}:")
        print(f"A({k}) = {a_k}")
        print(f"f({k}) = A({k}) - f({k-1}) = {a_k} - {f_prev} = {current_f}")
        
        f_prev = current_f

    print(f"\nThe total number of such sets T of size m = {m} is {final_answer}.")
    
    # Required final answer format
    print(f"<<<{final_answer}>>>")


if __name__ == '__main__':
    # You can change the values of n and m here
    # n must be a positive integer
    # m must be a non-negative integer
    n = 4
    m = 4
    solve_and_print(n, m)