import math

def solve():
    """
    Solves for the smallest positive integer M.
    
    The problem is equivalent to finding the maximum minimal range W of partial sums,
    such that W <= 2 * (1 - 10**(-M)).
    
    The worst-case sequence that maximizes this minimal range is shown to be
    a_1 = 1, a_k = 1 - delta for k > 1. This gives a minimal range of 2 - 2*delta.
    A specific "hardest" sequence gives W = 2 * (1 - 1/2**(N-1)).
    
    So, we solve the inequality:
    2 * (1 - 1/2**(N-1)) <= 2 * (1 - 10**(-M))
    1 - 1/2**(N-1) <= 1 - 10**(-M)
    10**(-M) <= 1/2**(N-1)
    -M <= log10(1/2**(N-1))
    -M <= -(N-1) * log10(2)
    M >= (N-1) * log10(2)
    
    We need to find the smallest integer M satisfying this.
    """
    N = 100000
    
    # Calculate the lower bound for M
    lower_bound = (N - 1) * math.log10(2)
    
    # M must be the smallest integer greater than or equal to this value
    M = math.ceil(lower_bound)
    
    print(f"Let N = 100000.")
    print(f"The condition is equivalent to M >= (N-1) * log10(2).")
    print(f"M >= (100000 - 1) * log10(2)")
    print(f"M >= 99999 * {math.log10(2)}")
    print(f"M >= {lower_bound}")
    print(f"The smallest positive integer M satisfying this is {M}.")
    
solve()