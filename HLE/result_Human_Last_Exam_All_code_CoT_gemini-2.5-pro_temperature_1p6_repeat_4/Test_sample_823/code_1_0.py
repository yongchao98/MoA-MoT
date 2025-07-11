import networkx as nx

def demonstrate_counterexamples():
    """
    This function demonstrates why options A, C, and E are not necessarily true
    by using the family of grid graphs as a counterexample.
    Grid graphs have bounded degree (at most 4) and unbounded treewidth.
    """

    print("Analyzing properties of grid graphs as a counterexample family.")
    print("The class of n x n grid graphs has bounded degree (<= 4) and unbounded treewidth (treewidth of an n x n grid is n).\n")

    grid_sizes = [5, 10, 20]

    for n in grid_sizes:
        print(f"--- Properties for a {n}x{n} grid graph ---")
        G = nx.grid_2d_graph(n, n)

        # --- Analysis for Option A: Induced Cycles ---
        # Grid graphs are bipartite, so they have no odd-length cycles.
        is_bipartite = nx.is_bipartite(G)
        print(f"Analysis for A (Induced Cycles):")
        print(f"  Is the {n}x{n} grid bipartite? {is_bipartite}.")
        if is_bipartite:
            print("  Since it's bipartite, it has no odd cycles. It cannot have induced cycles for odd lengths k=3, 5, ...")
        
        # --- Analysis for Option C: Expander Graphs ---
        # The Cheeger constant of n x n grids tends to 0 as n -> infinity.
        # We calculate an upper bound for it by splitting the grid vertically.
        partition_size = n * (n // 2)
        cut_size = n
        cheeger_upper_bound = cut_size / partition_size
        print(f"\nAnalysis for C (Expander Graphs):")
        print(f"  A specific cut gives a Cheeger constant upper bound: h(G) <= |E(S,V-S)|/|S|")
        print(f"  For a vertical bisection, h(G) <= {cut_size} / {partition_size} = {cheeger_upper_bound:.4f}")
        print(f"  This value approaches 0 as n grows, so grids are not an expander family.")

        # --- Analysis for Option E: Clique-Minors ---
        # Grid graphs are planar, which means their largest clique-minor is bounded by 4.
        # networkx has a planarity check, though it can be slow for large graphs.
        # Note: nx.check_planarity might not be in older versions or might require a specific installation setup.
        # The fact that grids are planar is a fundamental property.
        print(f"\nAnalysis for E (Clique Minors):")
        print(f"  Is the {n}x{n} grid planar? True (by definition).")
        print("  Since it's planar, it cannot contain a K5 minor. The clique-minor size is bounded by 4.")
            
        print("-" * (32 + len(str(n))*2))

if __name__ == '__main__':
    try:
        demonstrate_counterexamples()
    except ImportError:
        print("Please install networkx library: pip install networkx")
