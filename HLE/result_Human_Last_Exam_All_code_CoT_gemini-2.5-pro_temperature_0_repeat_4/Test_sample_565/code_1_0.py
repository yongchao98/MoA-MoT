import networkx as nx
import collections

def solve():
    """
    Calculates the number of non-isomorphic vertex-transitive graphs on 8
    vertices for each degree j from 0 to 7.
    """

    # This list contains functions, each constructing one of the 12 known
    # non-isomorphic vertex-transitive graphs on 8 vertices.
    graph_constructors = [
        # n_0: The empty graph on 8 vertices (degree 0)
        lambda: nx.empty_graph(8),
        
        # n_1: A perfect matching, i.e., 4 disjoint edges (degree 1)
        lambda: nx.from_edgelist([(0, 1), (2, 3), (4, 5), (6, 7)]),
        
        # n_2: The cycle graph on 8 vertices (degree 2)
        lambda: nx.cycle_graph(8),
        # n_2: The disjoint union of two 4-cycles (degree 2)
        lambda: nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4)),
        
        # n_3: The cubical graph (graph of the cube) (degree 3)
        lambda: nx.cubical_graph(),
        # n_3: The disjoint union of two complete graphs K4 (degree 3)
        lambda: nx.disjoint_union(nx.complete_graph(4), nx.complete_graph(4)),
        
        # n_4: The complement of the cubical graph (degree 4)
        lambda: nx.complement(nx.cubical_graph()),
        # n_4: The complete bipartite graph K_{4,4} (degree 4)
        lambda: nx.complete_bipartite_graph(4, 4),
        
        # n_5: The complement of the 8-cycle (degree 5)
        lambda: nx.complement(nx.cycle_graph(8)),
        # n_5: The complement of the disjoint union of two 4-cycles (degree 5)
        lambda: nx.complement(nx.disjoint_union(nx.cycle_graph(4), nx.cycle_graph(4))),
        
        # n_6: The complement of the perfect matching (Cocktail party graph) (degree 6)
        lambda: nx.complement(nx.from_edgelist([(0, 1), (2, 3), (4, 5), (6, 7)])),
        
        # n_7: The complete graph K8 (degree 7)
        lambda: nx.complete_graph(8),
    ]

    # Use a dictionary to store the counts n_j for each degree j
    counts = collections.defaultdict(int)

    for constructor in graph_constructors:
        g = constructor()
        # For a regular graph, all vertices have the same degree.
        # We can get the degree from any vertex.
        # For an empty graph, the degree is 0.
        if g.number_of_nodes() > 0:
            any_node = next(iter(g.nodes()))
            degree = g.degree(any_node)
        else:
            degree = 0
        
        counts[degree] += 1

    # Prepare the final list [n_0, n_1, ..., n_7]
    result_list = [counts[j] for j in range(8)]

    # Print the final result in the specified format
    print(result_list)

solve()