import networkx as nx
from networkx.algorithms import isomorphism

def demonstrate_graph_properties():
    """
    Analyzes the properties of a class of graphs with bounded degree and unbounded treewidth
    using the family of n x n grid graphs as a concrete example.
    """
    # Let C be the class of n x n grid graphs. For any n, the max degree is 4 (a constant).
    # The treewidth of an n x n grid is n, so the class has unbounded treewidth.
    n = 10
    G = nx.grid_2d_graph(n, n) # Nodes are tuples (x, y) with 0 <= x,y < n
    print(f"--- Analyzing the {n}x{n} grid graph as an example ---")
    print(f"This graph class has bounded degree (<=4) and unbounded treewidth (tw = n).")

    # A. For each k, there is a graph in C containing an induced cycle of length k.
    print("\n--- Analysis of Option A: Induced Cycles ---")
    # Grid graphs are bipartite, meaning they only contain even-length cycles.
    # Therefore, they do not contain induced cycles of odd lengths like 3 or 5.
    print("Grids are bipartite, so they have no odd-length cycles (e.g., C3, C5).")
    print("Thus, a class of grids does not contain induced cycles for all lengths k.")
    print("Conclusion: Option A is FALSE.")

    # B. For each k, there is a graph in C containing the k-by-k-grid as a subgraph.
    print("\n--- Analysis of Option B: Grid Subgraphs ---")
    # The Grid Minor Theorem states that large treewidth implies a large grid minor.
    # For graphs with bounded degree, this implies a large grid subgraph.
    # Our example class (grids) trivially satisfies this.
    k = 4
    grid_k = nx.grid_2d_graph(k, k)
    matcher = isomorphism.GraphMatcher(G, grid_k)
    is_subgraph = matcher.subgraph_is_isomorphic()
    print(f"The {n}x{n} grid contains a {k}x{k} grid as a subgraph: {is_subgraph}")
    print("Theoretically, for any class C with bounded degree and unbounded treewidth, this must be true.")
    print("Conclusion: Option B MUST BE TRUE.")

    # C. C is a family of expander graphs.
    print("\n--- Analysis of Option C: Expander Graphs ---")
    # Expanders have high connectivity. Grids do not, they have sparse cuts.
    # We can cut the grid in half vertically.
    S = set((i, j) for i in range(n // 2) for j in range(n))
    cut_size = len(nx.edge_boundary(G, S))
    # The Cheeger constant is roughly cut_size / |S|. For expanders, this is bounded from below.
    # For an n x n grid, the cut is n edges for a set of size n*n/2.
    # The ratio is n / (n^2 / 2) = 2/n.
    print(f"A cut separating the first {n//2} rows from the rest has {cut_size} edges.")
    print(f"The ratio |cut|/|S| is {cut_size / len(S):.2f}, which approaches 0 as n increases.")
    print("Grids are not an expander family.")
    print("Conclusion: Option C is FALSE.")

    # D. For each k, there is a graph in C containing an induced matching of size k.
    print("\n--- Analysis of Option D: Induced Matchings ---")
    # A large grid subgraph (from B) implies a large induced matching can be constructed.
    # Example: select horizontal edges that are spaced far apart.
    matching = []
    size_goal = 5
    for i in range(n // 4):
        if len(matching) >= size_goal: break
        u, v = (4*i, 0), (4*i+1, 0)
        if G.has_edge(u, v):
            matching.append((u, v))
    # We confirmed theoretically that this construction is an induced matching.
    print(f"Found an induced matching of size {len(matching)}: {matching}")
    print("Since we can find arbitrarily large grids (Option B), we can find arbitrarily large induced matchings.")
    print("Conclusion: Option D is also TRUE, as a consequence of B.")

    # E. The graphs in C contain clique-minors of unbounded size.
    print("\n--- Analysis of Option E: Clique-Minors ---")
    # Grids are planar graphs. By Kuratowski's theorem, they cannot contain a K5 minor.
    # Therefore, the size of clique-minors in grids is bounded (by 4).
    K5 = nx.complete_graph(5)
    matcher_K5 = isomorphism.GraphMatcher(G, K5)
    print(f"Is K5 a subgraph of the {n}x{n} grid? {matcher_K5.subgraph_is_isomorphic()}")
    print("Grids are planar and cannot contain a K5-minor, so clique-minor size is bounded.")
    print("Conclusion: Option E is FALSE.")
    
    print("\n--- Final Verdict ---")
    print("Options B and D must both be true. However, the presence of an arbitrarily large grid subgraph (B)")
    print("is the fundamental structural property that characterizes graphs with unbounded treewidth and bounded degree.")
    print("The presence of a large induced matching (D) is a weaker property that follows from (B).")
    print("Therefore, B is the most precise and fundamental answer.")

if __name__ == '__main__':
    demonstrate_graph_properties()