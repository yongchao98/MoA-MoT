import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.

    Args:
        n: A positive integer.
        m: A positive integer.

    Returns:
        The minimum possible diameter.
    """
    # A tree with V > 2 vertices must have at least 2 leaves.
    # n+2 > 2 => n > 0. Since n is a positive integer, this is always true.
    # So we must have m >= 2, except for the case of 2 vertices (n=0, but n is positive).
    if n == 1 and m == 1:
        # Special case: n+2=3, m=1 is not a valid tree (must have 2 leaves)
        # However, if we allowed this, it would be a malformed graph.
        # Based on standard tree properties for V>=3, m must be >= 2.
        # We will assume n and m are such that a valid tree can be formed.
        pass

    if m > n:
        # This implies m = n + 1, as the number of internal nodes I = n-m+2 must be >= 1.
        # With I=1, the tree is a star graph with diameter 2.
        diameter = 2
    elif m == n:
        # I = n - m + 2 = 2. The core is two vertices connected by an edge.
        # The diameter is 3 (leaf-core-core-leaf).
        diameter = 3
    else: # m < n
        # I = n - m + 2 > 2. The core has more than 2 vertices.
        if 2 * m >= n + 1:
            # We have enough leaves to make the core tree a star graph (diameter 2).
            # The total diameter of G is 2 (for T_I) + 2 (for leaves) = 4.
            diameter = 4
        elif n == 2 * m:
            # The core structure leads to a diameter of 5.
            diameter = 5
        else: # n > 2 * m
            # The core structure must be more "string-like", leading to a diameter of 6.
            diameter = 6
            
    return diameter

def solve():
    """
    Solves the problem for a given n and m.
    You can change the values of n and m here to test different cases.
    """
    # Let's use an example where n and m are provided.
    # For instance, let's take n = 5, m = 3.
    n = 5
    m = 3

    # Explanation of the case for n=5, m=3:
    # V = n+2 = 7 vertices.
    # L = m = 3 leaves.
    # The condition is m < n (3 < 5).
    # We check 2*m vs n+1: 2*3 = 6, n+1 = 5+1 = 6.
    # So, 2*m >= n+1 holds (6 >= 6).
    # According to our case analysis, the minimum diameter should be 4.
    
    diameter = get_minimum_diameter(n, m)
    
    print(f"For n = {n} and m = {m}:")
    print(f"The number of vertices is n + 2 = {n+2}.")
    print(f"The number of leaves is m = {m}.")
    print(f"The minimum possible value for the diameter of G is {diameter}.")

solve()
