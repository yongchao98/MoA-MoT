def solve_rank():
    """
    Calculates the rank of the kernel of the Hurewicz homomorphism for the given space Y.
    Y is the connected sum of X1, X2, and X3.
    """
    # Orders of the cyclic groups pi_1(X_i)
    n1 = 5
    n2 = 8
    n3 = 2
    
    # The group G = pi_1(Y) is the free product Z_n1 * Z_n2 * Z_n3
    # The kernel K of the Hurewicz map is the commutator subgroup [G, G].
    # The index of K in G is the order of the abelianization of G.
    index_K = n1 * n2 * n3
    
    print(f"The fundamental group is G = Z_{n1} * Z_{n2} * Z_{n3}.")
    print("The kernel K is the commutator subgroup [G,G].")
    print(f"The index of K in G is |G/[G,G]| = |Z_{n1} x Z_{n2} x Z_{n3}| = {n1} * {n2} * {n3} = {index_K}.")
    
    # Using Bass-Serre theory, the rank of the free group K can be calculated.
    # rank(K) = 1 - |V(K\\T)| + |E(K\\T)|
    # T is the Bass-Serre tree of the graph of groups for G.
    # The graph of groups has a central vertex v0 with group {1}, and leaf vertices
    # v1, v2, v3 with groups G1, G2, G3.
    
    # Number of vertices in the quotient graph K\T
    # |V| = |K\G/{1}| + |K\G/G1| + |K\G/G2| + |K\G/G3|
    V_v0 = index_K
    V_v1 = index_K / n1
    V_v2 = index_K / n2
    V_v3 = index_K / n3
    num_vertices = V_v0 + V_v1 + V_v2 + V_v3
    
    print("\nCalculating the number of vertices in the quotient graph K\\T:")
    print(f"|V| = [G:K] + [G:KG1] + [G:KG2] + [G:KG3]")
    print(f"    = {index_K} + {index_K}/{n1} + {index_K}/{n2} + {index_K}/{n3}")
    print(f"    = {int(V_v0)} + {int(V_v1)} + {int(V_v2)} + {int(V_v3)} = {int(num_vertices)}")
    
    # Number of edges in the quotient graph K\T
    # The edge groups are all trivial. There are 3 edges in the graph of groups.
    # |E| = 3 * |K\G/{1}| = 3 * [G:K]
    num_edges = 3 * index_K
    
    print("\nCalculating the number of edges in the quotient graph K\\T:")
    print(f"|E| = 3 * [G:K] = 3 * {index_K} = {num_edges}")

    # Rank formula for the free group K
    rank = 1 - num_vertices + num_edges
    
    print("\nThe rank of K as a free group is given by rank(K) = 1 - |V| + |E|.")
    print(f"rank(K) = 1 - {int(num_vertices)} + {int(num_edges)} = {int(rank)}")

solve_rank()