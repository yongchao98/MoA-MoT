def solve_max_n(k):
    """
    Calculates the maximum value of n for a given k based on the formula n = 2k - 1.

    The problem asks for the maximum value of n such that there exists a
    k-uniform intersecting family F of subsets of {1, ..., n} with
    full differences of size k-1.

    The derivation shows:
    1. For n < 2k - 1, any two k-sets intersect in at least 2 elements,
       making differences of size k-1 impossible.
    2. For n = 2k - 1, the family of all k-sets F = C(n, k) is intersecting
       and can be shown to have full differences of size k-1.
    3. For n > 2k - 1, it can be shown that no such family exists.

    Therefore, the maximum value of n is 2k - 1.
    """
    if not isinstance(k, int) or k < 1:
        print("Error: k must be a positive integer.")
        return

    n = 2 * k - 1
    
    print(f"For k = {k}, the maximum value of n is given by the formula n = 2*k - 1.")
    print("Substituting the value of k:")
    print(f"n = 2*{k} - 1")
    print(f"n = {2*k} - 1")
    print(f"n = {n}")

# You can test the function with any integer value of k >= 1.
# For example, let's use k=3.
solve_max_n(3)