import sys

def solve():
    """
    This function calculates the minimum possible diameter of a tree G with n+2 vertices and m leaves.
    The values of n and m are read from the user.
    """
    # The problem is given in terms of variables n and m.
    # We will first define them, for example:
    n = 10
    m = 4

    # This is an example of user input.
    # In a real scenario, you might read this from stdin, e.g.
    # if len(sys.argv) > 2:
    #     n = int(sys.argv[1])
    #     m = int(sys.argv[2])
    # else:
    #    print("Please provide n and m as command-line arguments.")
    #    # Set default values for demonstration
    #    n = 10
    #    m = 4
    
    print(f"Given n = {n} and m = {m}:")
    print(f"The tree G has {n+2} vertices and {m} leaves.")
    
    # A tree must have at least 2 leaves if it has more than 2 vertices.
    # Also, the number of leaves m can be at most n+1 (for a star graph on n+2 vertices).
    if m < 2 or m > n + 1:
        if n+2 > 2 and m < 2:
            print("Invalid input: A tree with more than 2 vertices must have at least 2 leaves.")
        elif m > n + 1:
            print("Invalid input: A tree with n+2 vertices can have at most n+1 leaves.")
        return

    # Case 1: The tree is a star graph. It has 1 internal node.
    # Number of internal nodes I = (n+2) - m = 1  => m = n+1
    if m == n + 1:
        min_diameter = 2
        print(f"Condition 'm = n+1' is met ({m} = {n}+1).")
        print(f"The tree is a star graph, so the minimum diameter is {min_diameter}.")
    
    # Case 2: The internal nodes form a path of length 1 (an edge).
    # Number of internal nodes I = (n+2) - m = 2  => m = n
    elif m == n:
        min_diameter = 3
        print(f"Condition 'm = n' is met ({m} = {n}).")
        print(f"The tree's internal nodes form a single edge, minimum diameter is {min_diameter}.")
    
    # Case 3: More than 2 internal nodes (m <= n-1).
    # We split this into two sub-cases.
    elif m <= n - 1:
        # Subcase 3a: There are "enough" leaves to keep the internal subgraph very compact (like a star).
        # The condition for this is 2*m >= n+1.
        if 2 * m >= n + 1:
            min_diameter = 4
            print(f"Condition 'm <= n-1' ({m} <= {n}-1) and '2m >= n+1' ({2*m} >= {n}+1) is met.")
            print(f"The minimum diameter is {min_diameter}.")
        # Subcase 3b: There are "not enough" leaves, forcing the internal subgraph to be more path-like.
        # The condition for this is 2*m < n+1.
        else:
            min_diameter = n - 2 * m + 5
            print(f"Condition 'm <= n-1' ({m} <= {n}-1) and '2m < n+1' ({2*m} < {n}+1) is met.")
            print(f"The minimum diameter is n - 2m + 5 = {n} - 2*{m} + 5 = {min_diameter}.")

solve()
