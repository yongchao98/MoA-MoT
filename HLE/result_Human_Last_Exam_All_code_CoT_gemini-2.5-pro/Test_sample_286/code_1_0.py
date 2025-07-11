import math

def solve():
    """
    Solves for the smallest positive integer M based on the problem description.
    
    The problem can be reformulated as finding the smallest positive integer M such that:
    sup_{a_i} min_{s_i} (span of partial sums) <= 2 * (1 - 10^{-M})
    
    Let N = 100000.
    The value of the supremum of the minimum span, W_N, for a sequence of length N is a result from discrepancy theory.
    For large N, a known bound which represents the worst-case is W_N = 2 - 2/N.
    
    So we must satisfy the inequality:
    2 - 2/N <= 2 * (1 - 10^{-M})
    1 - 1/N <= 1 - 10^{-M}
    10^{-M} <= 1/N
    -M <= -log10(N)
    M >= log10(N)
    
    We need to find the smallest integer M satisfying this condition.
    """
    
    N = 100000
    
    # We need to find the smallest integer M such that M >= log10(N).
    # This is equivalent to ceil(log10(N)).
    log_N = math.log10(N)
    M = math.ceil(log_N)
    
    print(f"Let N be the number of terms, N = {N}.")
    print("The condition to be satisfied is M >= log10(N).")
    print(f"M >= log10({N})")
    print(f"M >= {log_N}")
    print(f"The smallest integer M satisfying this condition is {M}.")

solve()