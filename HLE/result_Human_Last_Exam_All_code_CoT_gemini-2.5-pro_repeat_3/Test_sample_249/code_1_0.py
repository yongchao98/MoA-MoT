import math

def solve():
    """
    Calculates the minimum possible diameter of a tree with n+2 vertices and m leaves.
    """
    # Read n and m from user input.
    # For demonstration, we will use some example values.
    # To use with your own values, uncomment the next two lines.
    # n = int(input("Enter n: "))
    # m = int(input("Enter m: "))
    
    # Example values:
    n = 5
    m = 3

    print(f"For n = {n} and m = {m}:")
    
    # The number of internal (non-leaf) vertices is I = (n+2) - m.
    internal_nodes = n + 2 - m

    if internal_nodes == 1:
        # This corresponds to m = n + 1.
        # The structure is a star graph, diameter is 2.
        diameter = 2
        print(f"The number of internal nodes is {internal_nodes}.")
        print("This corresponds to a star graph, so the minimum diameter is 2.")

    elif internal_nodes == 2:
        # This corresponds to m = n.
        # The structure is two internal nodes connected, with leaves attached to them.
        # The diameter is 3.
        diameter = 3
        print(f"The number of internal nodes is {internal_nodes}.")
        print("The internal nodes form a path of length 1, so the minimum diameter is 3.")

    elif internal_nodes >= 3:
        # This corresponds to m <= n - 1.
        print(f"The number of internal nodes is {internal_nodes}.")
        if 2 * m >= n + 1:
            # We have enough leaves to make the internal nodes a star-like structure.
            # The diameter is 4.
            diameter = 4
            print(f"Since 2*m >= n+1 ({2*m} >= {n+1}), a compact star-like structure for internal nodes is possible.")
            print("The minimum diameter is 4.")
        else:
            # We don't have enough leaves for a star-like structure.
            # This forces a more path-like structure for the internal nodes.
            # The diameter is n + 5 - 2m.
            diameter = n + 5 - 2 * m
            print(f"Since 2*m < n+1 ({2*m} < {n+1}), the internal nodes must form a more elongated structure.")
            print(f"The minimum diameter is calculated as n + 5 - 2*m.")
            print(f"Result: {n} + 5 - 2 * {m} = {diameter}")
    else:
        # This case (I <= 0) should not happen for a valid tree with positive n, m.
        # A tree with V vertices has at most V-1 leaves. m <= n+1. So I = n+2-m >= 1.
        print("Invalid input: The number of leaves is too high for the number of vertices.")
        return

solve()