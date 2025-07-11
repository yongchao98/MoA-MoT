import itertools

def solve_m31_codim2():
    """
    Calculates and explains the number of codimension 2 boundary strata of the moduli space M_bar_{3,1}.

    The codimension of a boundary stratum corresponds to the number of nodes on the curve,
    which is equal to the number of edges in its dual graph. We are looking for stable dual
    graphs with g=3, n=1, and k=2 edges.

    The genus of a stable curve is given by the formula:
    g = h^1(G) + sum(g_v for v in V(G))
    where h^1(G) = |E| - |V| + 1 is the number of loops in the graph.

    The stability condition for each vertex v (representing an irreducible component) is:
    2*g_v - 2 + n_v > 0
    where g_v is the genus of the component and n_v is the number of special points
    (nodes and marked points) attached to it.
    """
    g = 3
    n = 1
    codimension = 2
    num_edges = codimension

    total_strata_count = 0
    
    print("This script calculates the number of codimension 2 boundary strata of M_bar_{3,1}.")
    print("This is equivalent to counting stable dual graphs with g=3, n=1, and 2 edges.")
    print("We proceed by classifying graphs by their number of vertices |V|.")
    print("-" * 50)

    # Case 1: |V| = 1 vertex
    num_v_1 = 1
    h1_1 = num_edges - num_v_1 + 1  # 2 - 1 + 1 = 2 loops
    sum_gv_1 = g - h1_1             # 3 - 2 = 1
    
    # Topology: 1 vertex with 2 loops. The marked point must be on this vertex.
    # Each loop creates one node.
    # n_v = 1 (marked point) + 2 (nodes from 2 loops) = 3
    g_v1 = sum_gv_1
    count_1 = 0
    if 2 * g_v1 - 2 + 3 > 0:
        count_1 = 1
    total_strata_count += count_1
    print(f"Case |V|=1: Graph has {h1_1} loops. Sum of vertex genera = {g} - {h1_1} = {sum_gv_1}.")
    print(f"  - Configuration: 1 vertex of genus {g_v1} with the marked point and 2 self-loops.")
    print(f"    - Stability check: 2*({g_v1}) - 2 + 3 = {2 * g_v1 - 2 + 3} > 0. This is a valid stratum.")
    print(f"  - Strata found: {count_1}")
    print("-" * 50)

    # Case 2: |V| = 2 vertices
    num_v_2 = 2
    h1_2 = num_edges - num_v_2 + 1 # 2 - 2 + 1 = 1 loop
    sum_gv_2 = g - h1_2            # 3 - 1 = 2
    
    # The graph must have one loop. Two topologies are possible.
    
    # Topology 2A: Loop on one vertex (v1), bridge to the other (v2).
    # Vertices v1 and v2 are distinguishable.
    strata_2a = set()
    for g_v1 in range(sum_gv_2 + 1):
        g_v2 = sum_gv_2 - g_v1
        # Leg on v1 (which has the loop): n1=3 (leg+loop+bridge), n2=1 (bridge)
        if (2*g_v1 - 2 + 3 > 0) and (2*g_v2 - 2 + 1 > 0):
            # Canonical representation: (genus of looped vtx, genus of other vtx, leg on looped vtx)
            strata_2a.add((g_v1, g_v2, True))
        # Leg on v2: n1=2 (loop+bridge), n2=2 (leg+bridge)
        if (2*g_v1 - 2 + 2 > 0) and (2*g_v2 - 2 + 2 > 0):
            strata_2a.add((g_v1, g_v2, False))
    count_2a = len(strata_2a)

    # Topology 2B: Two edges connecting v1 and v2.
    # Vertices are indistinguishable if their genera are equal.
    strata_2b = set()
    # Assume leg is on v1 by symmetry: n1=3 (leg+2 bridges), n2=2 (2 bridges)
    for g_v1 in range(sum_gv_2 + 1):
        g_v2 = sum_gv_2 - g_v1
        if (2*g_v1 - 2 + 3 > 0) and (2*g_v2 - 2 + 2 > 0):
            # Canonical representation: sorted tuple of (genus, leg_present)
            strata_2b.add(tuple(sorted([(g_v1, True), (g_v2, False)])))
    count_2b = len(strata_2b)
    
    count_2 = count_2a + count_2b
    total_strata_count += count_2
    print(f"Case |V|=2: Graph has {h1_2} loop. Sum of vertex genera = {g} - {h1_2} = {sum_gv_2}.")
    print(f"  - Topology A (loop on one vertex, bridge between): {count_2a} strata.")
    print(f"  - Topology B (two vertices connected by two edges): {count_2b} strata.")
    print(f"  - Strata found: {count_2}")
    print("-" * 50)

    # Case 3: |V| = 3 vertices
    num_v_3 = 3
    h1_3 = num_edges - num_v_3 + 1 # 2 - 3 + 1 = 0 loops (tree)
    sum_gv_3 = g - h1_3            # 3 - 0 = 3
    
    # Only topology is a chain: v1 -- v2 -- v3
    strata_3 = set()
    
    # Partitions of 3 into 3 parts
    genus_partitions = {p for i in range(sum_gv_3 + 1) for j in range(sum_gv_3 - i + 1) for p in set(itertools.permutations((i, j, sum_gv_3 - i - j)))}

    # Subcase 3a: Leg on an end vertex (e.g., v1).
    # n_v1=2, n_v2=2, n_v3=1
    count_3a = 0
    stable_configs_3a = set()
    for g_v1, g_v2, g_v3 in genus_partitions:
        if (2*g_v1 - 2 + 2 > 0) and (2*g_v2 - 2 + 2 > 0) and (2*g_v3 - 2 + 1 > 0):
            # Canonical representation: (middle genus, sorted tuple of end genera)
            stable_configs_3a.add((g_v2, tuple(sorted((g_v1, g_v3)))))
    count_3a = len(stable_configs_3a)

    # Subcase 3b: Leg on the middle vertex (v2).
    # n_v1=1, n_v2=3, n_v3=1
    count_3b = 0
    stable_configs_3b = set()
    for g_v1, g_v2, g_v3 in genus_partitions:
        if (2*g_v1 - 2 + 1 > 0) and (2*g_v2 - 2 + 3 > 0) and (2*g_v3 - 2 + 1 > 0):
             stable_configs_3b.add((g_v2, tuple(sorted((g_v1, g_v3)))))
    count_3b = len(stable_configs_3b)

    count_3 = count_3a + count_3b
    total_strata_count += count_3
    print(f"Case |V|=3: Graph is a tree (h1=0). Sum of vertex genera = {g} - {h1_3} = {sum_gv_3}.")
    print(f"  - Topology is a chain of 3 vertices.")
    print(f"  - Strata with leg on an end vertex: {count_3a}.")
    print(f"  - Strata with leg on the middle vertex: {count_3b}.")
    print(f"  - Strata found: {count_3}")
    print("-" * 50)

    print("Final Calculation:")
    print(f"Total number of strata = (from |V|=1) + (from |V|=2) + (from |V|=3)")
    print(f"Total number of strata = {count_1} + {count_2} + {count_3} = {total_strata_count}")
    return total_strata_count

if __name__ == '__main__':
    final_answer = solve_m31_codim2()
    print(f"<<<{final_answer}>>>")
