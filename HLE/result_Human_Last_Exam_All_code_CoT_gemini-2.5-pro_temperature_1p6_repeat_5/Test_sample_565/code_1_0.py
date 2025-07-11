import networkx as nx

def solve():
    """
    This function calculates the number of isomorphism classes of vertex-transitive
    graphs with 8 vertices for each degree from 0 to 7.
    """

    # This list will store one representative from each isomorphism class.
    # The graphs are identified based on established graph theory results.
    vt_graphs = []

    # Degree 0: The empty graph on 8 vertices.
    vt_graphs.append(nx.empty_graph(8))

    # Degree 1: A perfect matching on 8 vertices (4 disjoint copies of K2).
    vt_graphs.append(nx.from_edgelist([(0, 1), (2, 3), (4, 5), (6, 7)]))

    # Degree 2: The 8-cycle and two disjoint 4-cycles.
    vt_graphs.append(nx.cycle_graph(8))
    vt_graphs.append(nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4)))

    # Degree 3: The cube graph, two circulants, and two disjoint 4-cliques.
    vt_graphs.append(nx.hypercube_graph(3))
    vt_graphs.append(nx.circulant_graph(8, [1, 4]))
    vt_graphs.append(nx.circulant_graph(8, [2, 4]))
    vt_graphs.append(nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4)))
    
    # Graphs for degrees 4, 5, 6, 7 are complements of graphs for degrees 3, 2, 1, 0.
    # We find all graphs up to degree 3, then take their complements to find the rest.
    graphs_for_complements = list(vt_graphs)
    graphs_for_complements.reverse() # To get complements for degrees 4..7 from graphs of degrees 3..0
    
    for g in graphs_for_complements:
        vt_graphs.append(nx.complement(g))

    # Initialize the count for each degree j.
    n_j = [0] * 8
    
    for g in vt_graphs:
        # A vertex-transitive graph is regular, so all nodes have the same degree.
        # We can find the degree of the graph from its first node.
        # This handles the empty graph case as well since it has nodes.
        degree = g.degree[0]
        if 0 <= degree < 8:
            n_j[degree] += 1
            
    print(n_j)

solve()