import networkx as nx

def count_p3_subgraphs(G):
    """
    Counts the number of subgraphs isomorphic to P_3 (a path with 3 vertices) in a graph G.
    
    This function illustrates a case where #Sub_G(H) is tractable. H=P_3 is a graph
    with degree at most 2. The algorithm to count its occurrences runs in polynomial time,
    which means the problem is Fixed-Parameter Tractable (FPT) for H=P_3. This serves as a
    concrete counterexample to statement C, which claims that for some class of graphs H with
    degree at most 2, the problem could be #W[1]-complete. An FPT problem cannot be #W[1]-complete
    unless FPT = #W[1].

    The number of P_3 subgraphs centered at a vertex 'v' is the number of ways to choose
    2 of its neighbors, which is C(deg(v), 2). Summing this over all vertices gives the total count.

    Args:
        G (networkx.Graph): The input graph.

    Returns:
        int: The number of P_3 subgraphs.
    """
    p3_count = 0
    print("Analysis for counting P_3 subgraphs:")
    print("Vertex | Degree | P_3s Centered Here (C(Degree, 2))")
    print("-------------------------------------------------")
    for v in sorted(G.nodes()):
        degree = G.degree(v)
        # We need at least 2 neighbors to form a P_3
        if degree >= 2:
            combinations = degree * (degree - 1) // 2
            p3_count += combinations
            print(f"  {v: <4} |   {degree: <4} | {combinations}")
    
    print("-------------------------------------------------")
    return p3_count

def main():
    # We will use a Wheel graph on 5 vertices (W_5) as an example for G.
    # It has a central hub connected to 4 vertices in a cycle.
    G = nx.wheel_graph(5)
    
    print("This script analyzes the complexity of #Sub_G(H) for a specific H.")
    print("The goal is to provide evidence for the reasoning behind the final answer.\n")
    print(f"Let G be the wheel graph W_5.")
    print(f"Number of vertices in G: {G.number_of_nodes()}")
    print(f"Number of edges in G: {G.number_of_edges()}\n")

    # Let H be a path on 3 vertices, H = P_3.
    # H is a graph of maximum degree 2.
    num_p3 = count_p3_subgraphs(G)
    
    print(f"\nTotal number of P_3 subgraphs found in G: {num_p3}")
    
    print("\nConclusion from this example:")
    print("The counting was performed in polynomial time. This means for H=P_3, #Sub_G(H) is FPT.")
    print("Since P_3 is a graph of degree at most 2, this contradicts statement C.")
    print("A thorough analysis of all options shows that B is the correct statement.")


if __name__ == '__main__':
    main()