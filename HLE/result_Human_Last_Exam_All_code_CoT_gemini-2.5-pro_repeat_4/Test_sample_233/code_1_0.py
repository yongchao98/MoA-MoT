def solve_topology_problem():
    """
    Calculates the smallest positive integer g for the topology problem.

    Let g_Sigma be the genus of the initial surface Sigma.
    Let g_C be the genus of the capping surface C.
    The genus of the final closed surface Sigma' is g = g_Sigma + g_C.

    The problem asks for the smallest g that works for ANY embedding of Sigma.
    This means we must find the g for the "worst-case" embedding.
    The worst-case is when the capping surface C must have the highest possible
    genus. This occurs when all of the handles of Sigma are topologically
    linked with its boundary component.

    The number of handles on a surface is equal to its genus.
    Each handle of Sigma that is linked with the boundary forces the genus of
    the capping surface C to increase by 1.
    """

    # The genus of the given surface Sigma.
    genus_Sigma = 10

    # The number of handles on Sigma is equal to its genus.
    num_handles_Sigma = genus_Sigma

    # In the worst-case embedding, all handles are linked with the boundary.
    # The required genus of the capping surface C is equal to the number of
    # linked handles it must avoid intersecting.
    max_genus_C = num_handles_Sigma

    # The genus g of the final surface Sigma' is the sum of the genus of
    # Sigma and the genus of the capping surface C for the worst case.
    g = genus_Sigma + max_genus_C

    print("The genus of the final surface, g, is the sum of the genus of the initial surface (Sigma) and the genus of the capping surface (C).")
    print(f"Genus of Sigma: {genus_Sigma}")
    print("The genus of the capping surface depends on the embedding of Sigma.")
    print(f"In the worst case, all {num_handles_Sigma} handles of Sigma are linked with its boundary.")
    print(f"This requires a capping surface of genus {max_genus_C} to avoid intersections.")
    print("\nThe final equation is:")
    print(f"g = genus_Sigma + max_genus_C")
    print(f"g = {genus_Sigma} + {max_genus_C}")
    print(f"g = {g}")
    print("\nTherefore, the smallest positive integer g is:")
    print(g)

solve_topology_problem()