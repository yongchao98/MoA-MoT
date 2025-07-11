import itertools

def count_k_cliques(graph, k):
    """
    Counts the number of k-cliques in a graph using a brute-force approach.
    A graph is represented by an adjacency list (a dictionary where keys are nodes
    and values are sets of their neighbors).

    This algorithm's runtime is approximately O(n^k), where n is the number of vertices.
    This is not fixed-parameter tractable (FPT). Statement B asserts that for the
    graph class G specified in the problem, no FPT algorithm is expected to exist
    because the problem is #W[1]-complete.

    Args:
        graph (dict): The graph as an adjacency list.
        k (int): The size of the clique to count.

    Returns:
        int: The number of k-cliques in the graph.
    """
    nodes = list(graph.keys())
    if k > len(nodes):
        return 0
    
    clique_count = 0
    # Iterate through all combinations of k vertices
    for potential_clique in itertools.combinations(nodes, k):
        is_clique = True
        # Check if an edge exists between every pair of vertices in the combination
        for u, v in itertools.combinations(potential_clique, 2):
            if u not in graph or v not in graph[u]:
                is_clique = False
                break
        if is_clique:
            clique_count += 1
            
    return clique_count

def main():
    """
    Main function to demonstrate clique counting and explain the answer.
    """
    # Let's use a complete graph on 5 vertices (K5) as an example.
    # The class of all graphs that are subgraphs of some clique is
    # the class of all graphs, which is somewhere dense and closed under subgraphs.
    # While the problem states G is not the class of all graphs, this example
    # serves to illustrate the #k-Clique problem itself.
    G = {
        0: {1, 2, 3, 4},
        1: {0, 2, 3, 4},
        2: {0, 1, 3, 4},
        3: {0, 1, 2, 4},
        4: {0, 1, 2, 3},
    }
    
    # We want to count subgraphs isomorphic to H, where H is a clique.
    # Let's count 3-cliques (triangles).
    k = 3
    
    num_cliques_found = count_k_cliques(G, k)
    
    # The number of k-cliques in a K_n graph is given by the binomial coefficient "n choose k".
    # For K5 and k=3, this is "5 choose 3" = (5*4*3)/(3*2*1) = 10.
    
    print("Analyzing Statement B: If H is the class of all cliques, then #Sub_G(H) is #W[1]-complete.")
    print("This statement concerns the problem of counting k-cliques in a graph G from the class G.")
    print(f"Let's demonstrate with an example: counting {k}-cliques in a 5-clique graph.")
    print(f"The number of {k}-cliques found by the algorithm is: {num_cliques_found}.")
    print("The theoretical value is '5 choose 3', which is 10.")
    print("\nBased on established theorems in parameterized complexity, the problem of counting k-cliques")
    print("remains #W[1]-complete on any somewhere dense class of graphs that is closed under subgraphs.")
    print("The properties of the class G in the problem statement match these conditions exactly.")
    print("Therefore, statement B is the correct choice.")

if __name__ == '__main__':
    main()
