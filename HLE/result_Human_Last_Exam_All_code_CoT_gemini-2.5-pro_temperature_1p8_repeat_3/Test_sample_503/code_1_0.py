import math

def solve():
    """
    Calculates the number of internal adjunctions from [23] to [37]
    in the simplex category Delta.
    """
    m = 23
    n = 37

    # The number of internal adjunctions is the number of order-preserving maps
    # f: [m] -> [n] such that f(0) = 0.
    # This is determined by the sequence of values f(1), ..., f(m) which must be
    # non-decreasing and take values in {0, ..., n}.
    #
    # This is a combination with repetition problem: choosing k items from a set
    # of size N with replacement.
    # Number of items to choose: k = m = 23
    # Number of options for each item: N = n + 1 = 38
    # The formula is C(N + k - 1, k)
    k = m
    N = n + 1
    
    # Calculate the binomial coefficient C(N + k - 1, k) = C(38 + 23 - 1, 23) = C(60, 23)
    num_adjunctions = math.comb(N + k - 1, k)
    
    print(f"The number of internal adjunctions from [{m}] to [{n}] is the number of order-preserving maps f such that f(0)=0.")
    print(f"This is equivalent to choosing {m} values from the {n+1} integers {{0, ..., {n}}} with replacement.")
    print(f"The number of ways is given by the binomial coefficient C({N}+{k}-1, {k}).")
    print(f"The final equation is C({N+k-1}, {k}).")
    print(f"C({N+k-1}, {k}) = {num_adjunctions}")

solve()