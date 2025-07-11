def solve_surface_genus():
    """
    Calculates the smallest positive integer g for the problem.

    Let Σ be a smoothly embedded oriented surface of genus 10 in R^3 with a single
    unknotted boundary component. We want to find the smallest positive integer g
    such that, regardless of our choice of Σ, there exists a smoothly embedded
    oriented closed surface Σ' of genus g with Σ ⊆ Σ'.

    The solution involves two main parts:
    1. The genus of the original surface, g(Σ).
    2. The genus of the "capping" surface, g(C), needed to make Σ a closed surface.

    The final genus g is the sum g(Σ) + g(C), maximized over all possible embeddings of Σ.
    """

    # The genus of the given surface Σ is 10.
    g_sigma = 10

    # The problem asks for a g that works for *any* embedding of Σ. We must consider
    # the "worst-case" embedding. A surface of genus 10 can be thought of as a disk
    # with 10 handles. These handles can be topologically linked with the unknotted
    # boundary of Σ.
    #
    # To form a closed surface Σ' = Σ ∪ C, the capping surface C must not intersect Σ.
    # Each handle of Σ that is linked with the boundary forces the capping surface C
    # to have a handle of its own to avoid the intersection.
    #
    # Since Σ has 10 handles, it's possible to construct an embedding where all 10
    # handles are linked with the boundary. This forces the minimal capping surface C
    # to have a genus of 10. The genus of Σ itself limits this complexity, so the
    # maximum required genus for the cap is 10.
    g_cap_max = 10

    # The genus of the resulting closed surface Σ' is the sum of the genera of the
    # two surfaces it is composed of. For the worst-case scenario, this is:
    # g = g(Σ) + g_cap_max
    g_final = g_sigma + g_cap_max

    print(f"The genus of the original surface Σ is {g_sigma}.")
    print(f"The maximum required genus for the capping surface C is {g_cap_max}.")
    print("The smallest integer g that works for any Σ is the sum for the worst-case scenario.")
    print(f"g = {g_sigma} + {g_cap_max} = {g_final}")

solve_surface_genus()