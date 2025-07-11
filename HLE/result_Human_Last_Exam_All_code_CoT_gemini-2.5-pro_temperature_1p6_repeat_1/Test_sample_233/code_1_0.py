def solve_genus_problem():
    """
    Calculates the smallest integer g for the topology problem.
    """
    # The genus of the initial surface Σ.
    g_sigma = 10

    # The problem asks for the smallest integer g such that for any choice of Σ,
    # there exists a closed surface Σ' of genus g containing Σ.
    #
    # Let g_u be the number of handles on one side of a dividing plane and g_d be
    # the number of handles on the other side. g_u + g_d = g_sigma.
    #
    # The minimum genus of a cap C required to close Σ is min(g_u, g_d).
    # The minimum genus of the resulting surface Σ' is g_sigma + min(g_u, g_d).
    #
    # We must find a g that works for all cases. This means g must be at least
    # the maximum of these minimum possible genera.
    
    # We want to find max(g_sigma + min(g_u, g_d)) over all g_u from 0 to g_sigma.
    # This is equivalent to g_sigma + max(min(g_u, g_sigma - g_u)).

    max_min_cap_genus = 0
    g_u_for_max = -1

    for g_u in range(g_sigma + 1):
        g_d = g_sigma - g_u
        min_cap_genus = min(g_u, g_d)
        if min_cap_genus > max_min_cap_genus:
            max_min_cap_genus = min_cap_genus
            g_u_for_max = g_u
            
    g_d_for_max = g_sigma - g_u_for_max

    # The final answer is the original genus plus this maximized minimum cap genus.
    final_g = g_sigma + max_min_cap_genus
    
    print(f"The initial surface Σ has a genus of {g_sigma}.")
    print("The 10 handles of Σ can be distributed on two sides of a plane.")
    print(f"Let g_u and g_d be the number of handles on each side, with g_u + g_d = {g_sigma}.")
    print("To form a closed surface, a cap must be added. The minimal genus of this cap is min(g_u, g_d).")
    print(f"The minimal genus of the resulting closed surface is {g_sigma} + min(g_u, g_d).")
    print("To find a genus g that works for any Σ, we must consider the worst-case scenario,")
    print("which is the one that requires the largest possible minimal genus.")
    print("This occurs when min(g_u, g_d) is maximized.")
    print(f"For a total of {g_sigma} handles, this happens when the handles are split as evenly as possible:")
    print(f"g_u = {g_u_for_max}, g_d = {g_d_for_max}")
    print(f"The maximum value for the minimal cap genus is min({g_u_for_max}, {g_d_for_max}) = {max_min_cap_genus}.")
    print("The final equation for the smallest possible integer g is:")
    print(f"{g_sigma} + {max_min_cap_genus} = {final_g}")

solve_genus_problem()