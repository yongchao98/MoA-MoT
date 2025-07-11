import networkx as nx

def solve_vertex_transitive_graph_count():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs with 8 
    vertices for each possible degree (0 to 7).

    The approach is constructive:
    1. Identify the non-isomorphic vertex-transitive graphs for degrees j = 0, 1, 2, 3
       based on known graph theory results.
    2. Construct these graphs using networkx.
    3. The graphs for degrees j = 4, 5, 6, 7 are the complements of the graphs
       for j = 3, 2, 1, 0, respectively. Construct these by taking complements.
    4. Count the number of graphs constructed for each degree.
    """
    
    # A dictionary to hold the lists of graphs for each degree
    graphs_by_degree = {j: [] for j in range(8)}

    # Degree 0: The empty graph on 8 vertices (E_8)
    g0 = nx.empty_graph(8)
    graphs_by_degree[0].append(g0)

    # Degree 1: Four disjoint edges (4K_2)
    g1 = nx.Graph()
    g1.add_nodes_from(range(8))
    g1.add_edges_from([(0, 1), (2, 3), (4, 5), (6, 7)])
    graphs_by_degree[1].append(g1)

    # Degree 2: The 8-cycle (C_8) and two disjoint 4-cycles (2C_4)
    g2_1 = nx.cycle_graph(8)
    g2_2 = nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))
    graphs_by_degree[2].extend([g2_1, g2_2])

    # Degree 3: The cubical graph (Q_3) and two disjoint K_4s (2K_4)
    g3_1 = nx.cubical_graph()
    g3_2 = nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4))
    graphs_by_degree[3].extend([g3_1, g3_2])

    # Degrees 4, 5, 6, 7: Complements of graphs for degrees 3, 2, 1, 0
    for j in range(4, 8):
        complement_degree = 7 - j
        for g in graphs_by_degree[complement_degree]:
            # The complement of a vertex-transitive graph is vertex-transitive
            graphs_by_degree[j].append(nx.complement(g))

    # Calculate the final list of counts n_j
    result = [len(graphs_by_degree[j]) for j in range(8)]

    # The problem asks to output the numbers in the final equation.
    # We will print the list in the specified format.
    print(result)

solve_vertex_transitive_graph_count()