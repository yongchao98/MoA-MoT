import math

def solve():
    """
    Calculates the number of internal adjunctions in the simplex category Delta
    from [n] to [m].
    """
    n = 23
    m = 37
    
    # The number of such adjunctions corresponds to the number of monotone maps
    # f: [n] -> [m] such that f(0) = 0.
    # This is equivalent to counting non-decreasing sequences of length n,
    # f(1), ..., f(n), where each element is in {0, ..., m}.
    # The formula for this is C(m+n, n).
    
    total = n + m
    
    # Calculate the binomial coefficient C(n+m, n)
    result = math.comb(total, n)
    
    print(f"The number of internal adjunctions from [{n}] to [{m}] is given by the binomial coefficient C({m}+{n}, {n}).")
    print(f"C({total}, {n}) = {result}")

solve()