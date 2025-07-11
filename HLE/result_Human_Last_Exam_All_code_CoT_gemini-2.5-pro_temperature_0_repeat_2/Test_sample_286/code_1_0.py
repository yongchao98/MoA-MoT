import math

def solve():
    """
    This function solves for the smallest positive integer M based on the problem description.
    
    Let N = 100000.
    The condition that such a sequence x_i exists is equivalent to the condition that
    the minimal possible "swing" (or range) of the partial sums is less than or equal to the
    length of the interval where the x_i must lie.
    
    Let S_k = sum_{i=1 to k} s_i * a_i, where s_i are signs (+1 or -1).
    The range of the sequence x_i is max(S_k) - min(S_k).
    The length of the allowed interval for x_i is (1 - 10^{-M}) - (-1 + 10^{-M}) = 2 - 2 * 10^{-M}.
    
    So, for any sequence a_i, we must be able to choose signs s_i such that:
    max(S_k) - min(S_k) <= 2 - 2 * 10^{-M}
    
    This means that the worst-case minimal range must satisfy the inequality:
    sup_{a} min_{s} (max(S_k) - min(S_k)) <= 2 - 2 * 10^{-M}
    
    This value, let's call it K_N, is known from discrepancy theory. For N even,
    K_N = 2 - 2/N.
    
    Here N = 100000.
    So, K_100000 = 2 - 2/100000 = 2 - 2e-5.
    
    The inequality is:
    2 - 2e-5 <= 2 - 2 * 10^{-M}
    -2e-5 <= -2 * 10^{-M}
    1e-5 >= 10^{-M}
    -5 >= -M
    5 <= M
    
    The smallest positive integer M is 5.
    """
    N = 100000
    
    # The supremum of the minimal range of partial sums for a sequence of length N is K_N.
    # For N even, K_N = 2 - 2/N.
    K_N = 2 - 2/N
    
    # We need K_N <= 2 - 2 * 10^{-M}.
    # 2 - 2/N <= 2 - 2 * 10^{-M}
    # -2/N <= -2 * 10^{-M}
    # 1/N >= 10^{-M}
    # log10(1/N) >= -M
    # -log10(N) >= -M
    # log10(N) <= M
    
    M = math.ceil(math.log10(N))
    
    print(f"The number of terms is N = {N}.")
    print(f"The condition is sup(min(range)) <= 2 - 2 * 10^(-M).")
    print(f"The value of sup(min(range)) for N={N} is K_N = 2 - 2/N.")
    print(f"So, we must have 2 - 2/{N} <= 2 - 2 * 10^(-M).")
    print(f"This simplifies to 1/{N} >= 10^(-M).")
    print(f"Taking log10 on both sides: -log10({N}) >= -M.")
    print(f"This is equivalent to M >= log10({N}).")
    print(f"log10({N}) = {math.log10(N)}.")
    print(f"Since M must be an integer, M is the ceiling of {math.log10(N)}.")
    print(f"The smallest integer M is {M}.")

solve()