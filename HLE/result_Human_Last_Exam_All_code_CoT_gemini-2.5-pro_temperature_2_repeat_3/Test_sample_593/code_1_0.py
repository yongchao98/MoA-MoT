import networkx as nx

def calculate_treewidth(graph):
    """Calculates the treewidth of a networkx graph."""
    # For some simple graph families, the treewidth is known.
    # networkx's generic algorithm can be slow, so we'll use known formulas.
    if nx.is_isomorphic(graph, nx.complete_graph(len(graph.nodes))):
        # Treewidth of K_n is n-1
        return len(graph.nodes) - 1
    if nx.is_isomorphic(graph, nx.cycle_graph(len(graph.nodes))):
        # Treewidth of C_n is 2 for n>=3
        return 2 if len(graph.nodes) >= 3 else len(graph.nodes) - 1
    # For other graphs, we could use the approximation or exact algorithm.
    # For this example, manual calculation is sufficient.
    # For G = C5_complement = C5, this is covered.
    return -1 # Fallback

def solve():
    """
    This function demonstrates the tightness of the treewidth bound for a graph union.
    It uses the example of a K5 graph decomposed into two C5 cycles.
    """
    
    # Let's consider an example to show the bound is tight.
    # We will construct H and G such that their union F=H U G is the complete graph K_5.
    # The treewidth of K_5 is 5 - 1 = 4.

    # Let H be a 5-cycle (C_5) on 5 vertices.
    t_H = 2  # The treewidth of a cycle C_n (n>=3) is 2.

    # Let G be the complement of C_5. The complement of C_5 is also a C_5.
    t_G = 2  # The treewidth of G is also 2.

    # H and G are on the same set of 5 vertices.
    # So, the intersection V(H) intersect V(G) is the set of all 5 vertices.
    k = 5

    # The union F = H U G is the complete graph K_5, since E(H) U E(G) = E(K_5).
    # The treewidth of F = K_5 is 5 - 1 = 4.
    tw_F = 4

    # The formula for the upper bound is max(t_H, t_G, k - 1).
    bound = max(t_H, t_G, k - 1)
    
    print("Demonstration of the tight upper bound for treewidth:")
    print("-----------------------------------------------------")
    print(f"Let H be a cycle on 5 vertices (C5). Its treewidth is t_H = {t_H}.")
    print(f"Let G be the complement of C5 (which is also a C5). Its treewidth is t_G = {t_G}.")
    print(f"The graphs share all their vertices, so the intersection size is k = {k}.")
    print(f"The resulting graph F = H U G is the complete graph K5, and its treewidth is {tw_F}.")
    
    print("\nThe tight upper bound is calculated as: max(t_H, t_G, k - 1).")
    print(f"The final equation for the bound is: max({t_H}, {t_G}, {k} - 1) = {bound}")
    
    print("\nSince the treewidth of F ({}) is equal to the calculated bound ({}), the bound is tight.".format(tw_F, bound))

solve()
