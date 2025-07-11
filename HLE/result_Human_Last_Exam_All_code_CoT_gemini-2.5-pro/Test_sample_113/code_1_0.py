def count_strata():
    """
    Calculates and explains the number of codimension 2 boundary strata of M_bar(3,1).

    This is equivalent to counting the number of stable decorated dual graphs with
    g=3, n=1, and k=2 nodes (edges). The script systematically enumerates these graphs
    based on the number of vertices (V).
    """
    total_strata = 0
    print("Finding the number of codimension 2 boundary strata of M_bar_{3,1}.\n")
    print("This corresponds to counting stable curves with 2 nodes.\n")
    print("The dual graphs of these curves must satisfy Σg_v = V, where V is the number of vertices.\n")

    # --- Case 1: V = 1 vertex ---
    print("--- Case: 1 Vertex (V=1) ---")
    # For V=1, the 2 edges must be loops on the single vertex.
    # Σg_v = V  =>  g_1 = 1
    # The vertex has 1 marked point (leg) and 2 loops (2*2=4 half-edges).
    # Degree d_1 = 1 + 4 = 5.
    # Stability: 2*g_1 - 2 + d_1 = 2*1 - 2 + 5 = 5 > 0. This is stable.
    count_v1 = 1
    print("Graph: A single component of genus 1 with two self-nodes (a genus 1 curve with two pairs of points identified) and the marked point.")
    print(f"  - Genus assignment (g_1) = (1)")
    print(f"  - Stability check for g=1, d=5: 2*1 - 2 + 5 = 5 > 0. OK.")
    print(f"Result for V=1: {count_v1} stratum\n")
    total_strata += count_v1

    # --- Case 2: V = 2 vertices ---
    print("--- Case: 2 Vertices (V=2) ---")
    # Σg_v = V  =>  g_1 + g_2 = 2
    # Two possible graph topologies: parallel edges or a chain with a loop.
    count_v2 = 0
    
    # Topology A: Two vertices connected by two parallel edges.
    print("Topology A: Two components connected at two points.")
    # Leg is on v1. Degrees: d_1=3, d_2=2.
    # Stability: 2*g_1 - 2 + 3 > 0 => g_1 >= 0.
    #            2*g_2 - 2 + 2 > 0 => g_2 >= 1.
    # Genus partitions (g1, g2) summing to 2 with g2>=1:
    print("  - Genus (1,1): A genus 1 curve attached to another genus 1 curve at two points, with the marked point on one. Stable. (1 stratum)")
    count_v2 += 1
    print("  - Genus (0,2): A genus 0 curve attached to a genus 2 curve at two points, with the marked point on the genus 0 curve. Stable. (1 stratum)")
    count_v2 += 1
    
    # Topology B: A chain with a loop (loop on v1, edge between v1 and v2).
    print("\nTopology B: An irreducible nodal curve attached to another curve.")
    # Subcase B1: Marked point on v1 (with the loop).
    # Degrees: d_1=4, d_2=1. Stability: g_1>=0, g_2>=1.
    print("  Subcase B1: Marked point on the component with the self-node.")
    print("    - Genus (1,1): A g=1 curve with a node, attached to another g=1 curve. Marked point on the first. Stable. (1 stratum)")
    count_v2 += 1
    print("    - Genus (0,2): A g=0 curve with a node, attached to a g=2 curve. Marked point on the first. Stable. (1 stratum)")
    count_v2 += 1

    # Subcase B2: Marked point on v2 (without the loop).
    # Degrees: d_1=3, d_2=2. Stability: g_1>=0, g_2>=1.
    print("  Subcase B2: Marked point on the component without the self-node.")
    print("    - Genus (1,1): A g=1 curve with a node, attached to another g=1 curve. Marked point on the second. Stable. (1 stratum)")
    count_v2 += 1
    print("    - Genus (0,2): A g=0 curve with a node, attached to a g=2 curve. Marked point on the second. Stable. (1 stratum)")
    count_v2 += 1
    
    print(f"Result for V=2: {count_v2} strata\n")
    total_strata += count_v2
    
    # --- Case 3: V = 3 vertices ---
    print("--- Case: 3 Vertices (V=3) ---")
    # Σg_v = V => g_1 + g_2 + g_3 = 3
    # The only connected graph with 3 vertices and 2 edges is a chain: v1 -- v2 -- v3.
    count_v3 = 0
    print("Topology: A chain of three components.")
    
    # Subcase C1: Marked point on an end vertex (v1).
    # Degrees: d_1=2, d_2=2, d_3=1.
    # Stability requires: g_1>=1, g_2>=1, g_3>=1.
    # The only genus partition is (1,1,1).
    print("  Subcase C1: Marked point on an end component.")
    print("    - Genus (1,1,1): A chain of three g=1 curves, marked point on one end. Stable. (1 stratum)")
    count_v3 += 1
    
    # Subcase C2: Marked point on the middle vertex (v2).
    # Degrees: d_1=1, d_2=3, d_3=1.
    # Stability requires: g_1>=1, g_3>=1, g_2>=0.
    print("  Subcase C2: Marked point on the central component.")
    # Genus partitions (g1,g2,g3) summing to 3 with g1,g3 >= 1:
    print("    - Genus (1,1,1): A chain of three g=1 curves, marked point in the middle. Stable. (1 stratum)")
    count_v3 += 1
    print("    - Genus {1,0,2}: A g=1 and g=2 curve attached to a central g=0 curve with the marked point. Stable. (1 stratum)")
    # Note: (1,0,2) and (2,0,1) for (g1,g2,g3) are isomorphic as decorated graphs.
    count_v3 += 1

    print(f"Result for V=3: {count_v3} strata\n")
    total_strata += count_v3
    
    # --- Final Calculation ---
    print("--- Total ---")
    print(f"Number of strata from V=1: {count_v1}")
    print(f"Number of strata from V=2: {count_v2}")
    print(f"Number of strata from V=3: {count_v3}")
    print(f"Total number of codimension 2 strata = {count_v1} + {count_v2} + {count_v3} = {total_strata}")

if __name__ == '__main__':
    count_strata()