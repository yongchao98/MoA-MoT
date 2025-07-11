def solve_surface_genus():
    """
    Calculates the smallest possible genus g for a closed surface Sigma'
    that contains a given surface Sigma.

    The problem provides the following information:
    - Sigma is a smoothly embedded oriented surface in R^3.
    - The genus of Sigma is 10.
    - Sigma has a single unknotted boundary component.
    - Sigma' is a smoothly embedded oriented closed surface.
    - Sigma is a subsurface of Sigma' (Sigma is contained in Sigma').
    """

    # Genus of the initial surface, Sigma.
    genus_sigma = 10

    # The condition Sigma_prime contains Sigma means that Sigma_prime can be formed
    # by attaching a "capping" surface, S_cap, to the boundary of Sigma.
    # The genus of Sigma_prime is the sum of the genera of its parts.
    # Equation: g(Sigma') = g(Sigma) + g(S_cap)

    # A key theorem in 3-manifold topology states that a boundary component
    # that is an unknot in R^3 can always be bounded by an embedded disk in the
    # complement of the surface. A disk is a surface of genus 0.
    # This means the minimum possible genus for the capping surface is 0,
    # regardless of the complexity of Sigma's embedding.
    min_genus_s_cap = 0

    # We calculate the required genus, g, using the minimal-genus capping surface.
    # This gives the minimal possible genus for Sigma_prime for any given Sigma.
    # g = genus_sigma + min_genus_s_cap
    
    # Since any surface containing Sigma must have a genus of at least 10,
    # and we have shown that a genus of 10 is always sufficient,
    # the smallest integer g that works for all such surfaces is 10.

    print("The genus of the initial surface Σ is g(Σ) = 10.")
    print("The final closed surface Σ' is constructed by capping the boundary of Σ.")
    print("The genus of Σ' is given by the equation: g(Σ') = g(Σ) + g(S_cap), where S_cap is the capping surface.")
    print("Because the boundary of Σ is unknotted, a fundamental theorem guarantees that we can always find a capping surface S_cap with genus 0.")
    print("Therefore, the required genus g is calculated as:")
    
    # Output the final equation with the numbers plugged in.
    print(f"g = {genus_sigma} + {min_genus_s_cap}")
    
    result = genus_sigma + min_genus_s_cap
    print(f"g = {result}")

solve_surface_genus()