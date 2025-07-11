def solve_moduli_problem():
    """
    Calculates the number of codimension 2 boundary strata of the moduli space M_bar_{3,1}.

    This is equivalent to counting the number of stable dual graphs with
    g=3, n=1, and k=2 edges.
    """
    g = 3
    n = 1
    num_edges = 2
    total_strata = 0

    print("We are counting the number of stable dual graphs for (g,n)=(3,1) with 2 edges.")
    print("A graph is stable if for every vertex v, 2*g_v - 2 + n_v > 0.")
    print("The total genus is given by g = h^1(G) + sum(g_v).\n")

    # --- Case 1: |V|=1 vertex, |E|=2 edges ---
    print("--- Case 1: Graph with 1 vertex and 2 edges (two loops) ---")
    num_vertices = 1
    h1 = num_edges - num_vertices + 1  # h^1(G) = 2 - 1 + 1 = 2
    # Genus formula: 3 = 2 + g_1  => g_1 = 1
    g1 = g - h1
    # Stability check: n_1 = 1 leg + 2*2 half-edges = 5
    n1 = n + 2 * num_edges
    stability_val = 2 * g1 - 2 + n1
    if stability_val > 0:
        count_case1 = 1
        print(f"Configuration: g_1={g1}. Stability: 2*{g1}-2+{n1}={stability_val} > 0. Valid.")
        print(f"Number of strata found: {count_case1}\n")
        total_strata += count_case1
    else:
        count_case1 = 0

    # --- Case 2: |V|=2 vertices, |E|=2 edges ---
    print("--- Case 2: Graph with 2 vertices and 2 parallel edges ---")
    num_vertices = 2
    h1 = num_edges - num_vertices + 1  # h^1(G) = 2 - 2 + 1 = 1
    # Genus formula: 3 = 1 + g_1 + g_2 => g_1 + g_2 = 2
    genus_sum = g - h1
    print(f"Constraint: g_1 + g_2 = {genus_sum}.")
    print("We check partitions of 2 into (g_1, g_2) satisfying stability.")
    # Point on v1: n1=3, n2=2. Stability: 2g1-2+3>0 (g1>=0), 2g2-2+2>0 (g2>=1)
    count_case2 = 0
    for g2_val in range(genus_sum + 1):
        if g2_val >= 1:  # Stability for v2
            g1_val = genus_sum - g2_val
            if g1_val >= 0:  # Stability for v1
                count_case2 += 1
                print(f"  - Valid partition: (g_1, g_2) = ({g1_val}, {g2_val}) with point on g_1 component.")
    print(f"Number of strata found: {count_case2}\n")
    total_strata += count_case2

    # --- Case 3: |V|=3 vertices, |E|=2 edges ---
    print("--- Case 3: Graph with 3 vertices in a path (v1-v2-v3) ---")
    num_vertices = 3
    h1 = num_edges - num_vertices + 1  # h^1(G) = 2 - 3 + 1 = 0
    # Genus formula: 3 = 0 + g_1 + g_2 + g_3 => g_1 + g_2 + g_3 = 3
    genus_sum = g - h1
    print(f"Constraint: g_1 + g_2 + g_3 = {genus_sum}.")

    # Subcase 3a: Point on an end vertex (v1)
    # n1=2, n2=2, n3=1. Stability: g1>=1, g2>=1, g3>=1
    # Only solution for sum=3 is (1,1,1)
    count_case3a = 1
    print(f"Subcase 3a (point on end vertex): requires g_1,g_2,g_3 >= 1.")
    print("  - Valid partition: (1, 1, 1).")
    print(f"Number of strata found: {count_case3a}")
    total_strata += count_case3a

    # Subcase 3b: Point on the middle vertex (v2)
    # n1=1, n2=3, n3=1. Stability: g1>=1, g2>=0, g3>=1
    count_case3b = 0
    print(f"Subcase 3b (point on middle vertex): requires g_1>=1, g_2>=0, g_3>=1.")
    # Iterate through g1, g3 >= 1
    for g1_val in range(1, genus_sum):
        for g3_val in range(1, genus_sum):
            g2_val = genus_sum - g1_val - g3_val
            if g2_val >= 0:
                # Avoid double counting symmetric cases like (1,0,2) and (2,0,1)
                if g1_val <= g3_val:
                    count_case3b += 1
                    print(f"  - Valid partition: ({g1_val}, {g2_val}, {g3_val}).")
    print(f"Number of strata found: {count_case3b}\n")
    total_strata += count_case3b

    # --- Final Summation ---
    print("--- Total ---")
    count_case3 = count_case3a + count_case3b
    print(f"Total number of strata = (from Case 1) + (from Case 2) + (from Case 3)")
    print(f"Total number of strata = {count_case1} + {count_case2} + {count_case3}")
    print(f"The final answer is {total_strata}.")

solve_moduli_problem()