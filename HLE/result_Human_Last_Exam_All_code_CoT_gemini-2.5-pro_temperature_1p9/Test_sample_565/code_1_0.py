import networkx as nx

def get_vtx_transitive_graphs_on_8_vertices():
    """
    This function generates all 14 non-isomorphic vertex-transitive graphs on 8 vertices.
    It does so by defining a base set of 7 graphs and generating their complements.
    This list is known from graph theory literature.
    """
    base_graphs = []
    
    # j=0: 8 isolated vertices (8K1)
    base_graphs.append(nx.empty_graph(8))
    
    # j=1: 4 disjoint edges (4K2)
    g1 = nx.Graph()
    g1.add_edges_from([(0, 1), (2, 3), (4, 5), (6, 7)])
    base_graphs.append(g1)

    # j=2: C8 and 2C4
    base_graphs.append(nx.cycle_graph(8))
    base_graphs.append(nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4)))
    
    # j=3: Q3, C8(1,4), and 2K4
    base_graphs.append(nx.cubical_graph())
    base_graphs.append(nx.circulant_graph(8, [1, 4]))
    base_graphs.append(nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4)))
    
    # A self-complementary graph on 8 vertices must have 14 edges.
    # For a j-regular graph, this means 8*j/2 = 14 => 4j=14, which has no integer solution for j.
    # Therefore, no graph on this list is self-complementary, and the full list of 14 graphs
    # is the union of the base set and their complements.
    
    all_graphs = []
    all_graphs.extend(base_graphs)
    for g in base_graphs:
        all_graphs.append(nx.complement(g))
        
    return all_graphs

def calculate_degree_counts():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8 vertices
    for each degree j from 0 to 7.
    """
    graphs = get_vtx_transitive_graphs_on_8_vertices()
    
    counts = [0] * 8
    
    for g in graphs:
        # For a regular graph, all nodes have the same degree. We can pick any node.
        if len(g.nodes) == 0:
            degree = 0
        else:
            degree = g.degree[0]
            
        if 0 <= degree <= 7:
            counts[degree] += 1
            
    return counts

# Calculate and print the result
counts = calculate_degree_counts()
print(f"[{counts[0]}, {counts[1]}, {counts[2]}, {counts[3]}, {counts[4]}, {counts[5]}, {counts[6]}, {counts[7]}]")