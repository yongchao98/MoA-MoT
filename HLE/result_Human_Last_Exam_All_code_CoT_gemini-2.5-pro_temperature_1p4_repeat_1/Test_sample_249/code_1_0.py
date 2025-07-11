import math

def get_minimum_diameter(n, m):
    """
    Calculates the minimum possible diameter of an undirected tree with
    n+2 vertices and m leaves.

    Args:
        n: A positive integer.
        m: A positive integer.

    Returns:
        The minimum possible diameter.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # A tree must have at least 2 leaves (for more than 1 vertex).
    # The number of leaves m cannot exceed the number of vertices minus one.
    if m < 2 or m > n + 1:
        print(f"Error: For n={n}, m must be between 2 and n+1.")
        return

    # The number of internal nodes (degree >= 2) is k.
    k = n + 2 - m
    
    diameter = 0

    # Case 1: 2*m >= n+1. The internal subgraph can be made very compact (a star).
    if 2 * m >= n + 1:
        if k == 1: # m = n+1, G is a star graph
            diameter = 2
            print(f"For n={n} and m={m}, the tree has k=1 internal node.")
            print(f"The minimum diameter is {diameter}.")

        elif k == 2: # m = n, T_I is an edge
            diameter = 3
            print(f"For n={n} and m={m}, the tree has k=2 internal nodes.")
            print(f"The core of these two nodes is an edge of length 1.")
            print(f"The diameter of the full tree is D(T_I) + 2 = 1 + 2 = {diameter}.")
        else: # k >= 3
            diameter = 4
            print(f"For n={n} and m={m}, the tree has k={k} (>=3) internal nodes.")
            print(f"Since 2m >= n+1, we can form a star-like core, D(T_I) = 2.")
            print(f"The diameter of the full tree is D(T_I) + 2 = 2 + 2 = {diameter}.")

    # Case 2: 2*m < n+1. The internal subgraph cannot be a star.
    else:
        # Special subcase: If m=2, the tree must be a path graph P_{n+2}.
        if m == 2:
            diameter = n + 1
            print(f"For n={n} and m={m}, the tree is a path of {n+2} vertices.")
            print(f"The diameter is {diameter}.")

        # General subcase for 2 < m < (n+1)/2
        else:
            # We form a balanced T_I with k nodes and m leaves to minimize its diameter.
            # The diameter of this T_I is 2*ceil((k-2)/m)+1.
            d_i = 2 * math.ceil((k - 2) / m) + 1
            diameter = d_i + 2
            print(f"For n={n} and m={m}, the number of internal nodes is k={k}.")
            print("The number of leaves m is not enough to form a star-like core.")
            print("The core T_I is a balanced tree with k nodes and m leaves.")
            print(f"The diameter of this core is D(T_I) = 2*ceil(({k}-2)/{m})+1 = {d_i}.")
            print(f"The total minimum diameter is D(T_I) + 2 = {d_i} + 2 = {diameter}.")
            
    return diameter

# You can change these values to test different scenarios
n_val = 10
m_val = 4
get_minimum_diameter(n_val, m_val)