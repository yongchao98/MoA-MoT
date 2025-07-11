def find_required_vertices_for_induced_matching(max_degree_d, matching_size_k):
    """
    Calculates the minimum number of vertices required to guarantee an
    induced matching of a certain size in a graph with a given maximum degree.

    The logic is based on the inequality: |M| >= n / (2*d),
    where |M| is the size of the induced matching, n is the number of vertices,
    and d is the maximum degree.

    To guarantee a matching of size k, we need |M| >= k.
    So, we need n / (2*d) >= k, which means n >= k * 2 * d.
    We return the smallest integer n that satisfies this.
    """

    if not isinstance(max_degree_d, int) or max_degree_d <= 0:
        print("Maximum degree d must be a positive integer.")
        return
    if not isinstance(matching_size_k, int) or matching_size_k <= 0:
        print("Desired matching size k must be a positive integer.")
        return

    # From the proof, the size of a maximal induced matching |M| is at least n / (2*d).
    # To ensure |M| is at least k, we need n / (2*d) >= k.
    # Therefore, n must be at least k * 2 * d.
    min_n = matching_size_k * 2 * max_degree_d

    print(f"Given a class of graphs C with maximum degree at most d = {max_degree_d} and unbounded treewidth:")
    print(f"To guarantee an induced matching of size at least k = {matching_size_k}, we can prove this is always possible.")
    print("\nProof reasoning:")
    print("1. A graph with n vertices and maximum degree d has an induced matching of size |M|.")
    # In the final print, we show the numbers involved in the equation.
    print(f"2. It is guaranteed that |M| >= n / (2 * {max_degree_d}).")
    print(f"3. We want |M| >= {matching_size_k}, so we need n / (2 * {max_degree_d}) >= {matching_size_k}.")
    print(f"4. This simplifies to n >= {matching_size_k} * 2 * {max_degree_d}, so n >= {min_n}.")
    print("\nConclusion:")
    print(f"Since C has unbounded treewidth, it contains graphs with an arbitrarily large number of vertices.")
    print(f"Thus, for k={matching_size_k}, we can always find a graph in C with at least {min_n} vertices, which guarantees an induced matching of size at least {matching_size_k}.")


# Example: Let d = 10 and we want to find a graph with an induced matching of size k = 100.
example_d = 10
example_k = 100
find_required_vertices_for_induced_matching(example_d, example_k)