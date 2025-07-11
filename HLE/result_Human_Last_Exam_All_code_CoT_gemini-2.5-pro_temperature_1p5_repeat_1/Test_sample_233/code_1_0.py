def solve_topology_problem():
    """
    This function calculates the smallest possible genus for the final surface Sigma'.
    """
    # The genus of the initial surface Sigma is given as 10.
    # This represents the number of "handles" on the surface.
    genus_sigma = 10

    # To create a closed surface (a surface without a boundary), we must "cap" or "patch"
    # the single boundary component of Sigma.
    # The simplest possible surface that can serve as a cap is a disk.
    # A disk has 0 handles, so its genus is 0.
    genus_cap = 0

    # The genus of the final, closed surface Sigma' is the sum of the genus of the
    # original surface and the genus of the capping surface.
    final_genus = genus_sigma + genus_cap

    print("The final closed surface, denoted as Σ', is formed by taking the original surface Σ and capping its boundary.")
    print("The genus of the final surface is the sum of the genus of the original surface and the genus of the cap.")
    print("The equation is: g(Σ') = g(Σ) + g(cap)")
    print(f"The calculation is: g = {genus_sigma} + {genus_cap}")
    print(f"The resulting smallest possible genus g is: {final_genus}")

solve_topology_problem()