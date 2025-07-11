import math

def solve():
    """
    Solves for the smallest positive integer M based on the problem description.
    
    Let N = 100000.
    The problem is to find the smallest positive integer M such that for any sequence
    a_1, ..., a_N in [0,1], there exists a sequence x_0, ..., x_N in
    [-1 + 10^{-M}, 1 - 10^{-M}] with |x_{i-1} - x_i| = a_i.
    
    This is equivalent to finding the supremum of the minimal range of partial sums.
    Let S_k = sum_{i=1 to k} s_i * a_i, where s_i is +1 or -1.
    We need to find signs s_i and a starting point x_0 such that all x_k = x_0 + S_k
    are in the given interval.
    
    By choosing x_0 = -(max(S_k) + min(S_k))/2, the sequence x_k is centered,
    and its values lie in [-R/2, R/2], where R = max(S_k) - min(S_k).
    
    We need [-R/2, R/2] to be a subset of [-1 + 10^{-M}, 1 - 10^{-M}].
    This means R/2 <= 1 - 10^{-M}.
    
    This must hold for the best choice of signs (which minimizes R), for any sequence a_i.
    Let R_min(a) be the minimal range for a given sequence 'a'.
    Let R* = sup_a(R_min(a)).
    The condition is R*/2 <= 1 - 10^{-M}.
    
    A result from discrepancy theory states that for N vectors in [0,1],
    R* = 2 * (1 - 1/N).
    
    So, (2 * (1 - 1/N)) / 2 <= 1 - 10^{-M}
    1 - 1/N <= 1 - 10^{-M}
    1/N >= 10^{-M}
    10^M >= N
    M >= log10(N)
    
    Since M must be an integer, M = ceil(log10(N)).
    """
    N = 100000
    
    # The supremum of the minimal range R* is 2 * (1 - 1/N)
    R_star = 2 * (1 - 1/N)
    
    # We need R*/2 <= 1 - 10^{-M}
    # (2 * (1 - 1/N)) / 2 <= 1 - 10^{-M}
    # 1 - 1/N <= 1 - 10^{-M}
    # 10^{-M} <= 1/N
    # -M <= -log10(N)
    # M >= log10(N)
    
    M_float = math.log10(N)
    M = math.ceil(M_float)
    
    print(f"Let N be the number of values a_i, so N = {N}.")
    print("The condition on M can be derived from the supremum of the minimal range of the partial sums, R*.")
    print(f"A known result from discrepancy theory states that R* = 2 * (1 - 1/N) = 2 * (1 - 1/{N}) = {R_star}.")
    print("The condition is R* / 2 <= 1 - 10^(-M).")
    print(f"Substituting R*, we get: {R_star / 2} <= 1 - 10^(-M)")
    print(f"This simplifies to: 1 - 1/{N} <= 1 - 10^(-M)")
    print(f"Which further simplifies to: 10^M >= N")
    print(f"Taking the base-10 logarithm: M >= log10(N)")
    print(f"M >= log10({N})")
    print(f"M >= {M_float}")
    print(f"Since M must be the smallest positive integer, we take the ceiling.")
    print(f"M = ceil({M_float}) = {M}")
    
solve()