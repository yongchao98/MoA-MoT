def solve_topology_problem():
    """
    Calculates the smallest positive integer g based on topological principles.

    The problem asks for the smallest integer g such that any smoothly embedded
    oriented surface Σ of genus 10 with a single unknotted boundary component
    can be contained within a smoothly embedded oriented closed surface Σ' of genus g.

    The logic is as follows:
    1. The containing surface Σ' is formed by Σ and a "capping surface" S,
       such that Σ' = Σ ∪ S.
    2. The genus of the combined surface is the sum of the individual genera:
       g(Σ') = g(Σ) + g(S).
    3. We are given g(Σ) = 10.
    4. The problem asks for a g that works for ANY choice of Σ. This means we must
       consider the worst-case embedding of Σ, which is the one that forces the
       capping surface S to have the maximum possible genus.
    5. A theorem in topology states that the maximum possible genus of the capping
       surface, g(S)_max, is equal to the genus of the surface being capped, g(Σ).
    6. Therefore, the required genus g is the sum of the genus of Σ and the
       maximum possible genus of the capping surface S.
    """
    # Genus of the initial surface Σ
    g_sigma = 10

    # According to a theorem by Scharlemann, the capping genus g(S) is bounded by g(Σ).
    # This bound is known to be sharp, meaning a "worst-case" surface exists
    # for which the capping genus is exactly g(Σ).
    # Therefore, the maximum required genus for the capping surface S is g(Σ).
    max_g_s = g_sigma

    # The smallest integer g that works for all possible surfaces Σ must be
    # large enough for this worst-case scenario.
    g = g_sigma + max_g_s

    # Output the reasoning and the final equation.
    print("The genus of the original surface is g(Σ) = 10.")
    print("The genus of the final closed surface is g(Σ') = g(Σ) + g(S), where S is the capping surface.")
    print("To find the value 'g' that works for any Σ, we must consider the worst-case scenario.")
    print("The worst-case Σ forces the capping surface S to have the maximum possible genus.")
    print("Topological theorems state that this maximum possible genus for S is equal to the genus of Σ.")
    print("Therefore, the calculation is:")
    print(f"{g_sigma} + {max_g_s} = {g}")

solve_topology_problem()