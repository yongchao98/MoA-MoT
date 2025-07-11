import itertools

def count_k_cliques(graph, k):
    """
    Counts the number of k-cliques in a graph.

    The graph is represented as an adjacency list (dictionary in Python).
    A k-clique is a subset of k vertices where every two distinct vertices
    in the subset are adjacent.

    This implementation demonstrates the #k-Clique problem from statement B.
    Its runtime is roughly O(k^2 * n^k), where n is the number of vertices,
    which is not fixed-parameter tractable.

    Args:
        graph (dict): The graph as an adjacency list.
                      Keys are vertices (e.g., integers),
                      values are lists of adjacent vertices.
        k (int): The size of the clique to count.

    Returns:
        int: The number of k-cliques in the graph.
    """
    nodes = list(graph.keys())
    if k > len(nodes):
        return 0

    clique_count = 0
    # Iterate over all combinations of k vertices
    for subset in itertools.combinations(nodes, k):
        is_clique = True
        # Check if every pair of vertices in the subset is connected
        for u, v in itertools.combinations(subset, 2):
            # We assume an undirected graph, so we check both directions
            # for adjacency list representation.
            if v not in graph.get(u, []) and u not in graph.get(v, []):
                is_clique = False
                break
        if is_clique:
            clique_count += 1
            
    return clique_count

if __name__ == '__main__':
    # Example Graph G (a diamond graph with an extra edge 0-3)
    # It has two 3-cliques: {0, 1, 2} and {0, 2, 3}
    # and one 4-clique: {0, 1, 2, 3} if we add all edges.
    G = {
        0: [1, 2, 3],
        1: [0, 2],
        2: [0, 1, 3],
        3: [0, 2]
    }
    k = 3
    
    num_cliques = count_k_cliques(G, k)
    
    print(f"The input graph G has vertices: {list(G.keys())}")
    print(f"The input graph G has edges defined by the adjacency list.")
    print(f"The parameter k (clique size) is: {k}")
    print(f"The number of subgraphs of G isomorphic to a {k}-clique is: {num_cliques}")
