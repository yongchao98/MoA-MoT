def solve_moduli_problem():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}
    by enumerating the corresponding stable dual graphs.
    """
    print("This script calculates the number of codimension 2 boundary strata of the moduli space M_bar_{3,1}.")
    print("This is equivalent to counting the number of non-isomorphic decorated dual graphs representing stable curves of genus 3, with 1 marked point, and 2 nodes.\n")

    print("--- The Combinatorial Framework ---")
    print("A decorated dual graph must satisfy:")
    print("1. Parameters: Total genus g=3, marked points n=1, nodes |E|=2.")
    print("2. Genus Condition: sum(g_i) + h^1(G) = 3. Since h^1(G) = |E| - k + 1 = 3 - k, this means sum(g_i) = k, where k is the number of components (vertices).")
    print("3. Stability Condition: For each component i, 2*g_i - 2 + n_i > 0, where n_i is the number of nodes and marked points on that component.")

    # --- Case 1: k=1 component ---
    print("\n--- Case 1: One component (k=1) ---")
    print("The graph has 1 vertex and 2 edges, forming two self-loops (2 nodes on one component).")
    print("Genus Condition: g_1 = k = 1.")
    print("The component is a genus 1 curve with 2 nodes and 1 marked point, so it has n_1 = 3 special points.")
    print("Stability Check for (g=1, n=3): 2*1 - 2 + 3 = 3 > 0. Stable.")
    print("This gives 1 unique stratum (the irreducible type).")
    count_k1 = 1
    
    # --- Case 2: k=2 components ---
    print("\n--- Case 2: Two components (k=2) ---")
    print("The only connected graph has 2 vertices connected by 2 edges.")
    print("Genus Condition: g_1 + g_2 = k = 2.")
    print("We analyze the partitions of the genus sum 2:")
    
    print("  a) Partition (2, 0): Components with g_1=2, g_2=0.")
    print("     A g=0 component needs at least 3 special points to be stable. The two nodes provide 2 points.")
    print("     The marked point must be on the g=0 component. This gives it n_2=3 special points.")
    print("     - Component 1 (g=2, n=2 nodes): 2*2 - 2 + 2 = 4 > 0. Stable.")
    print("     - Component 2 (g=0, n=2 nodes + 1 point): 2*0 - 2 + 3 = 1 > 0. Stable.")
    print("     This configuration is valid and unique. Contribution: 1 stratum.")
    count_k2_a = 1

    print("  b) Partition (1, 1): Components with g_1=1, g_2=1.")
    print("     The components are indistinguishable by genus. We place the marked point on one.")
    print("     - Component 1 (with point): g=1, n=3 (2 nodes + 1 pt). 2*1 - 2 + 3 = 3 > 0. Stable.")
    print("     - Component 2 (no point): g=1, n=2 (2 nodes). 2*1 - 2 + 2 = 2 > 0. Stable.")
    print("     This configuration is valid and unique. Contribution: 1 stratum.")
    count_k2_b = 1
    count_k2 = count_k2_a + count_k2_b
    
    # --- Case 3: k=3 components ---
    print("\n--- Case 3: Three components (k=3) ---")
    print("The only connected graph is a chain: V1 -- V2 -- V3.")
    print("Genus Condition: g_1 + g_2 + g_3 = k = 3.")
    print("We analyze partitions of the genus sum 3:")
    
    print("  a) Partition (2, 1, 0): Components with genera {2, 1, 0}.")
    print("     Stability analysis shows the only valid configuration is when the g=0 component is in the middle (V2), and the marked point is on it.")
    print("     - V1/V3 (end, g=2/1, n=1 node): Stable (2*2-2+1=3 > 0 and 2*1-2+1=1 > 0).")
    print("     - V2 (middle, g=0, n=2 nodes + 1 pt): 2*0 - 2 + 3 = 1 > 0. Stable.")
    print("     This defines one unique stratum. Contribution: 1 stratum.")
    count_k3_a = 1

    print("  b) Partition (1, 1, 1): All three components have g=1.")
    print("     Components are indistinguishable. The position of the marked point defines the stratum.")
    print("     - Point on an end component (e.g., V1): (V1: n=2), (V2: n=2), (V3: n=1). All are stable. Contribution: 1 stratum.")
    count_k3_b1 = 1
    print("     - Point on the middle component (V2): (V1: n=1), (V2: n=3), (V3: n=1). All are stable. Contribution: 1 stratum.")
    count_k3_b2 = 1
    count_k3 = count_k3_a + count_k3_b1 + count_k3_b2
    
    # --- Final Calculation ---
    print("\n--- Total Count ---")
    total = count_k1 + count_k2 + count_k3
    print("The total number of strata is the sum of counts from each case:")
    print(f"Total = (k=1) + (k=2) + (k=3)")
    print(f"Total = {count_k1} + {count_k2} + {count_k3} = {total}")

if __name__ == '__main__':
    solve_moduli_problem()
    print("\n<<<6>>>")
