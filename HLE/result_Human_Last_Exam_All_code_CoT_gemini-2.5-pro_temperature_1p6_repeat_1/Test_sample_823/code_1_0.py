import networkx as nx

def analyze_grid_expansion(k):
    """
    Analyzes the expansion property of a k x k grid graph.

    It calculates the ratio of the boundary size to the set size for a
    specific cut (the top half of the grid). In an expander family, this
    ratio should be bounded below by a positive constant.
    """
    print(f"--- Analyzing {k}x{k} Grid ---")
    
    # A k x k grid has k*k vertices. Max degree is 4. Treewidth is k.
    # The class of grid graphs has bounded degree and unbounded treewidth.
    G = nx.grid_2d_graph(k, k)
    n_nodes = G.number_of_nodes()

    # Define a set S as the top half of the grid's vertices
    # Let S contain the first k/2 rows
    num_rows_in_S = k // 2
    if num_rows_in_S == 0:
        print("Grid is too small to analyze this cut.")
        return

    S = {(i, j) for i in range(num_rows_in_S) for j in range(k)}
    size_S = len(S)

    # The boundary E(S, V\S) consists of the k horizontal edges connecting
    # row (num_rows_in_S - 1) with row (num_rows_in_S).
    boundary_edges_count = k

    # The ratio is an upper bound on the Cheeger constant.
    ratio = boundary_edges_count / size_S
    
    # Equation: ratio = |E(S, V\S)| / |S|
    print(f"Number of vertices in the grid |V|: {n_nodes}")
    print(f"We choose a set S (the top {num_rows_in_S} rows).")
    print(f"Size of the set |S|: {size_S}")
    print(f"Number of edges on the boundary |E(S, V\\S)|: {boundary_edges_count}")
    print("The expansion ratio for this cut is given by the equation:")
    print(f"ratio = |E(S, V\\S)| / |S|")
    print(f"ratio = {boundary_edges_count} / {size_S} = {ratio:.4f}")
    if k > 10:
        print(f"As k increases, this ratio (proportional to {k}/{k*k/2} = 2/{k}) approaches 0.")


print("Demonstrating that the class of grid graphs is not an expander family (Counterexample for Option C).\n")
analyze_grid_expansion(10)
analyze_grid_expansion(50)
analyze_grid_expansion(100)
