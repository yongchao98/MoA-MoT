import math

def solve():
    """
    This problem asks for the largest number n for a specific decomposition of a metric continuum X.
    The properties of X are:
    1. It contains 5 distinct points S = {a, b, c, d, e}.
    2. No proper subcontinuum of X contains any three of these points.

    The properties of the decomposition X = A_1 U ... U A_n are:
    1. Each A_i is a subcontinuum.
    2. Each A_i has a non-empty part not covered by the other A_j's.

    From these properties, we can deduce that each A_i must be a proper subcontinuum.
    If A_i were equal to X, then for any other A_j, A_j would be a subset of X = A_i,
    which would make A_j \\ (U_{k!=j} A_k) empty, a contradiction.

    Since each A_i is a proper subcontinuum, property (2) of X implies that each A_i
    can contain at most 2 points from the set S.

    To maximize n, we should associate each continuum A_i with the smallest possible
    distinguishing feature related to S. The most natural choice is to have each A_i
    contain exactly one unique pair of points from S.

    The problem thus reduces to finding the maximum number of unique pairs we can choose
    from a set of 5 points. This is a combinatorial problem, "5 choose 2".

    The formula for combinations is C(n, k) = n! / (k! * (n-k)!).
    """
    
    n = 5  # Number of points
    k = 2  # Size of the groups (pairs)
    
    # Calculate the number of combinations
    result = math.comb(n, k)
    
    # The equation is C(5, 2) = 5! / (2! * (5-2)!)
    print("The problem is equivalent to finding the number of distinct pairs of points from a set of 5 points.")
    print("This is calculated using the combination formula C(n, k) = n! / (k! * (n-k)!).")
    print(f"Here, n = {n} and k = {k}.")
    print(f"So we calculate C({n}, {k}) = {n}! / ({k}! * ({n}-{k})!)")
    print(f"C(5, 2) = 5! / (2! * 3!) = (5 * 4) / (2 * 1) = 10")
    print(f"The largest number n is: {result}")

solve()