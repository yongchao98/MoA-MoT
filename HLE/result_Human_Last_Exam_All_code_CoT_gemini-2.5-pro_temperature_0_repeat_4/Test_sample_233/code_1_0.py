def solve_surface_genus():
    """
    Calculates the smallest positive integer g for the problem.
    """
    # The genus of the initial surface Sigma is given.
    g_sigma = 10
    print(f"The genus of the initial surface, g(Sigma), is {g_sigma}.")

    # The final surface Sigma' is formed by capping Sigma with a surface C.
    # The genus of the final surface is g(Sigma') = g(Sigma) + g(C).
    # We need to find the maximum necessary genus for the cap, g(C),
    # to handle the worst-case embedding of Sigma.

    # The number of handles on Sigma determines the maximum possible complexity
    # of its embedding with respect to its boundary. A surface of genus 10
    # has 10 handles.
    num_handles_sigma = 10

    # In the worst-case scenario, all of Sigma's handles are topologically
    # linked with its boundary. To create a non-intersecting cap, the cap C
    # must have one handle for each linked handle of Sigma it needs to avoid.
    g_cap_max = num_handles_sigma
    print(f"The maximum necessary genus for the capping surface, g(C), is {g_cap_max}.")

    # The smallest integer g that works for any choice of Sigma is the one
    # that can accommodate this worst-case scenario.
    # g = g(Sigma) + g(C)_max
    g_final = g_sigma + g_cap_max

    print(f"The final genus g is calculated as: g = {g_sigma} + {g_cap_max}")
    print(f"The result of the equation is: {g_final}")

solve_surface_genus()