import math

def solve():
    """
    Calculates the number of internal adjunctions in the simplex category Delta
    from [n] to [m].
    """
    n = 23
    m = 37

    # The number of internal adjunctions from [n] to [m] is given by the
    # number of order-preserving maps f: [n] -> [m] such that f(0) = 0.
    # This is a combinatorial problem of choosing n values for f(1),...,f(n)
    # from the set {0, 1, ..., m} with replacement.
    # The formula for this is the binomial coefficient C(m+n, n).

    # Calculate the binomial coefficient C(60, 23)
    result = math.comb(m + n, n)

    # Print the equation and the result
    print(f"The number of internal adjunctions is C({m} + {n}, {n}).")
    print(f"C({m+n}, {n}) = {result}")

solve()