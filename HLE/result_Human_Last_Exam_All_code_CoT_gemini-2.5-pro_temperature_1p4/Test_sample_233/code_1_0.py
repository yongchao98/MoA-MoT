def solve_genus_problem():
    """
    This function calculates the smallest positive integer g for the problem.
    The logic is based on the topological properties of the surfaces.

    g_sigma: The genus of the initial surface Σ.
    g_cap: The genus of the "capping" surface C.
    g_final: The genus of the final closed surface Σ'.

    The relationship is g_final = g_sigma + g_cap.
    We are given g_sigma = 10.
    We need to find the g_final that works for *any* such Σ, which means we must
    consider the worst-case Σ.
    A worst-case Σ is one that is topologically "linked" with its own boundary.
    This structure prevents capping with a simple disk (g_cap = 0).
    The obstruction can be overcome by using a cap with one handle (a punctured torus),
    which has g_cap = 1. This is sufficient for all cases.
    So, the maximum necessary g_cap is 1.
    """
    g_sigma = 10
    # Maximum required genus for the capping surface in the worst-case scenario.
    g_cap_max_required = 1

    g_final = g_sigma + g_cap_max_required

    # Output the equation and the final result.
    print(f"The genus of the original surface Σ is {g_sigma}.")
    print(f"The required genus for the capping surface C in the worst-case is {g_cap_max_required}.")
    print("The genus g of the final surface Σ' is g(Σ) + g(C).")
    print(f"{g_sigma} + {g_cap_max_required} = {g_final}")

solve_genus_problem()