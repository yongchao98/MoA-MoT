def sum_of_powers(n, p):
    """Calculates sum of p-th powers from 1 to n."""
    if p == 1:
        return n * (n + 1) // 2
    if p == 2:
        return n * (n + 1) * (2 * n + 1) // 6
    if p == 3:
        return (n * (n + 1) // 2)**2
    return -1

def b1_N_g(g):
    """Computes the first l2-betti number of the group N_g."""
    return 2 * g**2 - 2 * g

def sum_b1_G_vi(i):
    """
    Computes the sum of Betti numbers for the vertex group G_{v_i} = *_{g=2 to i} N_g.
    This corresponds to 2 * sum_{g=2 to i} (g^2 - g).
    """
    if i < 2:
        return 0
    # sum_{g=2 to i} (g^2 - g) = (sum_{g=1 to i} g^2 - sum_{g=1 to i} g) - (1^2 - 1)
    sum_g_sq = sum_of_powers(i, 2)
    sum_g = sum_of_powers(i, 1)
    return 2 * (sum_g_sq - sum_g)

def calculate_total_betti():
    """
    Calculates the first l2-betti number based on the plan.
    """
    # 1. Betti number of the graph Gamma
    b1_gamma = 16

    # 2. Sum of Betti numbers of vertex groups
    b1_gv1 = b1_N_g(100)
    
    sum_b1_other_vertices = 0
    for i in range(2, 16):
        sum_b1_other_vertices += sum_b1_G_vi(i)
        
    sum_b1_vertex_groups = b1_gv1 + sum_b1_other_vertices

    # 3. Sum of Betti numbers of edge groups
    # As per the plan, we assume all 30 edge groups are N_2.
    num_edges = 30
    b1_ge = b1_N_g(2)
    sum_b1_edge_groups = num_edges * b1_ge

    # 4. Final Calculation
    final_betti = b1_gamma + sum_b1_vertex_groups - sum_b1_edge_groups
    
    print("The first l2-Betti number is computed using the formula:")
    print("b1(Gamma) + sum(b1(G_v)) - sum(b1(G_e))")
    print("\nComponent values:")
    print(f"b1(Gamma) = {b1_gamma}")
    print(f"sum(b1(G_v)) = {sum_b1_vertex_groups}")
    print(f"sum(b1(G_e)) = {sum_b1_edge_groups}")
    print("\nFinal Equation:")
    print(f"{b1_gamma} + {sum_b1_vertex_groups} - {sum_b1_edge_groups} = {final_betti}")
    
    return final_betti

# Execute the calculation and store the final answer.
final_answer = calculate_total_betti()
# print(f"\nFinal Answer: {final_answer}")
# The problem asks for the answer in a specific format at the end.
# So we will prepare the final answer string.
# <<<29216>>>