def solve_surface_problem():
    """
    Calculates the smallest integer g for the topology problem.

    This function follows the logic that the final genus g must be achievable
    for all possible embeddings of the surface Sigma, including the "worst-case"
    embedding.
    """
    # The genus of the initial surface Σ.
    g_sigma = 10
    
    # The number of handles on a surface is equal to its genus.
    num_handles = g_sigma
    
    print("Let g_sigma be the genus of the initial surface Σ.")
    print(f"Given: g_sigma = {g_sigma}")
    print("\nA closed surface Σ' is formed by attaching a capping surface C.")
    print("The genus of the final surface Σ' is g' = g_sigma + g(C), where g(C) is the genus of the cap.")
    
    print("\nThe embedding of Σ can force the cap C to have a non-zero genus.")
    print("Specifically, handles of Σ can be linked with the boundary ∂Σ.")
    print("Let k be the number of handles linked with the boundary.")
    
    # The minimum required genus of the cap, g(C)_min, is k.
    # The number of linked handles k can range from 0 to the total number of handles.
    print(f"The number of handles is {num_handles}, so k can range from 0 to {num_handles}.")
    
    print("\nTo find a genus g' that works for ANY choice of Σ, we must consider the worst-case scenario.")
    print("The worst case is the surface that most restricts the choice of g', which occurs when k is maximum.")
    
    # In the worst case, all handles are linked with the boundary.
    k_worst = num_handles
    print(f"The maximal value for k is {k_worst}.")
    
    # For this worst-case surface, the minimum genus of the cap C is k_worst.
    g_c_min_worst_case = k_worst
    print(f"The minimum genus of the cap for this surface is g(C)_min = {g_c_min_worst_case}.")

    # Calculate the smallest possible genus for Σ' in this worst case.
    g_final_min = g_sigma + g_c_min_worst_case

    print("\nFor this worst-case surface, the smallest possible genus for Σ' is calculated as follows:")
    print(f"g' = g_sigma + g(C)_min")
    print(f"g' = {g_sigma} + {g_c_min_worst_case} = {g_final_min}")

    print(f"\nThus, for the worst-case Σ, the set of achievable genera for Σ' is {{{g_final_min}, {g_final_min+1}, ...}}.")
    print("For any other Σ (with k < 10), the set of achievable genera will start at a smaller number (10+k), but will also contain all integers from 20 onwards.")
    print("The intersection of all these sets is therefore {20, 21, 22, ...}.")
    print("The smallest integer g in this intersection is the answer.")
    print(f"\nSmallest positive integer g: {g_final_min}")

solve_surface_problem()