def solve_m31_codim2():
    """
    Calculates the number of codimension 2 boundary strata of M_bar_{3,1}.
    
    This corresponds to enumerating the topological types of stable curves
    of genus 3 with 1 marked point and 2 nodes.
    """
    strata_count = 0
    
    print("Enumerating codimension 2 boundary strata of M_bar(3,1):\n")

    # --- Case 1: 1 component (irreducible curve), V=1 ---
    # Graph: 1 vertex, 2 self-loops.
    # h^1 = |E| - |V| + 1 = 2 - 1 + 1 = 2.
    # g_total = g_comp + h^1 => 3 = g_comp + 2 => g_comp = 1.
    # k_v = 1 (marked point) + 2 (nodes) = 3.
    # Note: each self-loop corresponds to one node.
    g_v = 1
    k_v = 1 + 2
    if 2 * g_v - 2 + k_v > 0:
        strata_count += 1
        print(f"Stratum {strata_count}: Irreducible curve with two nodes")
        print("  - Graph: 1 vertex (component), 2 self-loops (nodes).")
        print(f"  - Component genus: {g_v}")
        print(f"  - Special points: 1 (marked point) + 2 (nodes) = {k_v}")
        print(f"  - Stability check: 2*g - 2 + k = 2*{g_v} - 2 + {k_v} = {2*g_v-2+k_v} > 0. Stable.")
        print("-" * 30)

    # --- Case 2: 2 components, V=2 ---
    # Graph: 2 vertices, 2 parallel edges.
    # h^1 = |E| - |V| + 1 = 2 - 2 + 1 = 1.
    # g_total = g1 + g2 + h^1 => 3 = g1 + g2 + 1 => g1 + g2 = 2.
    # Genus partitions (g1, g2): (1,1) and (2,0).
    
    # Subcase 2.1: g1=1, g2=1
    g = [1, 1]
    # Marked point on one component (symmetric case). k = [1(mp)+2(nodes), 2(nodes)] = [3,2]
    k = [3, 2]
    if (2 * g[0] - 2 + k[0] > 0) and (2 * g[1] - 2 + k[1] > 0):
        strata_count += 1
        print(f"Stratum {strata_count}: Two genus 1 components")
        print("  - Graph: 2 vertices, 2 parallel edges.")
        print("  - Component genera: (g1, g2) = (1, 1)")
        print(f"  - C1 (with m.p.): g=1, k=3. Stability: 2*1 - 2 + 3 = 3 > 0. Stable.")
        print(f"  - C2 (no m.p.):  g=1, k=2. Stability: 2*1 - 2 + 2 = 2 > 0. Stable.")
        print("-" * 30)
        
    # Subcase 2.2: g1=2, g2=0
    g = [2, 0]
    # Marked point on g=2 component: k=[3,2]. C2(g=0, k=2) is unstable (2*0-2+2=0).
    # Marked point on g=0 component: k=[2,3]
    k = [2, 3]
    if (2 * g[0] - 2 + k[0] > 0) and (2 * g[1] - 2 + k[1] > 0):
        strata_count += 1
        print(f"Stratum {strata_count}: A genus 2 and a genus 0 component")
        print("  - Graph: 2 vertices, 2 parallel edges.")
        print("  - Component genera: (g1, g2) = (2, 0), marked point on the genus 0 component.")
        print(f"  - C1 (no m.p.): g=2, k=2. Stability: 2*2 - 2 + 2 = 4 > 0. Stable.")
        print(f"  - C2 (with m.p.): g=0, k=3. Stability: 2*0 - 2 + 3 = 1 > 0. Stable.")
        print("-" * 30)
        
    # --- Case 3: 3 components, V=3 ---
    # Graph: A chain (path graph). v1 -- v2 -- v3.
    # h^1 = |E| - |V| + 1 = 2 - 3 + 1 = 0.
    # g_total = g1 + g2 + g3 + h^1 => g1 + g2 + g3 = 3.
    # Node connections: k_nodes = [1, 2, 1] for (v1, v2, v3).
    
    # Subcase 3.1: Genera (1,1,1)
    g = [1, 1, 1]
    # M.p. on an end component (v1): k=[1+1, 2, 1]=[2,2,1].
    k = [2, 2, 1]
    if (2*g[0]-2+k[0]>0) and (2*g[1]-2+k[1]>0) and (2*g[2]-2+k[2]>0):
        strata_count += 1
        print(f"Stratum {strata_count}: Chain of three genus 1 components, m.p. on an end")
        print("  - Graph: 3 vertices in a line.")
        print("  - Component genera: (g1, g2, g3) = (1, 1, 1). Marked point on an end component.")
        print(f"  - C1 (end, w/ mp): g=1, k=2. Stability: 2*1 - 2 + 2 = 2 > 0. Stable.")
        print(f"  - C2 (mid): g=1, k=2. Stability: 2*1 - 2 + 2 = 2 > 0. Stable.")
        print(f"  - C3 (end): g=1, k=1. Stability: 2*1 - 2 + 1 = 1 > 0. Stable.")
        print("-" * 30)
        
    # M.p. on the middle component (v2): k=[1, 2+1, 1]=[1,3,1].
    k = [1, 3, 1]
    if (2*g[0]-2+k[0]>0) and (2*g[1]-2+k[1]>0) and (2*g[2]-2+k[2]>0):
        strata_count += 1
        print(f"Stratum {strata_count}: Chain of three genus 1 components, m.p. in middle")
        print("  - Graph: 3 vertices in a line.")
        print("  - Component genera: (g1, g2, g3) = (1, 1, 1). Marked point on the middle component.")
        print(f"  - C1 (end): g=1, k=1. Stability: 2*1 - 2 + 1 = 1 > 0. Stable.")
        print(f"  - C2 (mid, w/ mp): g=1, k=3. Stability: 2*1 - 2 + 3 = 3 > 0. Stable.")
        print(f"  - C3 (end): g=1, k=1. Stability: 2*1 - 2 + 1 = 1 > 0. Stable.")
        print("-" * 30)

    # Subcase 3.2: Genera partition (2,1,0).
    # For stability with m.p. on middle comp, end comps must have g>0. So g_v=(1,0,2).
    # Marked point on middle (g=0) component: k=[1, 2+1, 1]=[1,3,1].
    g = [1, 0, 2]
    k = [1, 3, 1]
    if (2*g[0]-2+k[0]>0) and (2*g[1]-2+k[1]>0) and (2*g[2]-2+k[2]>0):
        strata_count += 1
        print(f"Stratum {strata_count}: Chain of components with genera 1, 0, 2")
        print("  - Graph: 3 vertices in a line.")
        print("  - Component genera: (g1, g2, g3) = (1, 0, 2). Marked point on the middle (g=0) component.")
        print(f"  - C1 (end): g=1, k=1. Stability: 2*1 - 2 + 1 = 1 > 0. Stable.")
        print(f"  - C2 (mid, w/ mp): g=0, k=3. Stability: 2*0 - 2 + 3 = 1 > 0. Stable.")
        print(f"  - C3 (end): g=2, k=1. Stability: 2*2 - 2 + 1 = 3 > 0. Stable.")
        print("-" * 30)
        
    print(f"\nThe total number of codimension 2 boundary strata is {strata_count}.")
    
    # We output the final calculation as a formula
    final_sum = "1 + 1 + 1 + 1 + 1 + 1 = 6"
    print("\nFinal calculation:")
    # Print each number in the equation.
    print(1, end="")
    for i in range(strata_count - 1):
        print(f" + {1}", end="")
    print(f" = {strata_count}")


solve_m31_codim2()