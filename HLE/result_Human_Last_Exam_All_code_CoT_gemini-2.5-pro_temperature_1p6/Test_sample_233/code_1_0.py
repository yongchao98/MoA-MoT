def solve_surface_genus():
    """
    Calculates the genus of the smallest containing closed surface.

    The problem asks for the smallest positive integer g such that for any smoothly
    embedded oriented surface Σ of genus 10 with a single unknotted boundary
    component in R^3, there exists a smoothly embedded oriented closed surface Σ'
    of genus g with Σ ⊆ Σ'.

    The genus of the final surface, g, is the sum of the genus of the original
    surface (g_sigma) and the genus of the "capping" surface (g_cap).
    g = g_sigma + g_cap

    The problem requires finding a g that works for any Σ. This means we must
    consider the worst-case embedding of Σ, which maximizes the necessary
    genus of the cap.

    A result from topology states that the maximum required genus for the
    capping surface is equal to the genus of the original surface.
    g_cap_max = g_sigma

    Therefore, the smallest g that is guaranteed to work for any Σ is:
    g = g_sigma + g_sigma = 2 * g_sigma
    """
    g_sigma = 10
    g_cap_max = g_sigma  # In the worst case, the cap needs the same genus.
    
    # The final genus g is the sum of the original surface's genus
    # and the worst-case cap's genus.
    g_final = g_sigma + g_cap_max
    
    print(f"The genus of the given surface Σ is g(Σ) = {g_sigma}.")
    print(f"The maximum necessary genus for the capping surface is g(S_fill) = {g_cap_max}.")
    print(f"The resulting genus of Σ' is g(Σ') = g(Σ) + g(S_fill) = {g_sigma} + {g_cap_max} = {g_final}.")

solve_surface_genus()
<<<20>>>