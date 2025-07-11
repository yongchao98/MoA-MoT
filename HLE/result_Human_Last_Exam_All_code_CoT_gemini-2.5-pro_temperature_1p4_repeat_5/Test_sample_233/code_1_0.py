import math

def solve_genus_problem():
    """
    Calculates the smallest possible genus for the final closed surface.

    The problem asks for the smallest positive integer g such that for any smoothly
    embedded oriented surface Sigma of genus 10 with a single unknotted boundary
    component, there exists a smoothly embedded oriented closed surface Sigma' of
    genus g such that Sigma is a subset of Sigma'.

    To form a closed surface Sigma' from a surface with a boundary Sigma, we must
    "cap" the boundary. The resulting genus is the sum of the genus of the original
    surface and the genus of the cap.

    g(Sigma') = g(Sigma) + g(cap)

    To find the smallest possible g(Sigma'), we must find the smallest possible g(cap).
    The genus of any surface is a non-negative integer. The smallest possible genus
    is therefore 0.

    The boundary of Sigma is an unknotted circle. A surface of genus 0 with one
    boundary component is a disk. By the Disk Theorem, an unknotted circle in R^3
    bounds a disk. Thus, we can always use a disk as the cap.
    """

    # The genus of the initial surface Sigma.
    genus_sigma = 10

    # The smallest possible genus for the capping surface (a disk) is 0.
    genus_cap = 0

    # The genus 'g' of the resulting closed surface Sigma' is the sum.
    g = genus_sigma + genus_cap

    print(f"The genus of the original surface Σ is {genus_sigma}.")
    print("To make the surface closed, we must cap its unknotted boundary.")
    print(f"The simplest cap is a disk, which has a genus of {genus_cap}.")
    print("The smallest possible genus g for the final surface Σ' is the sum of these genera.")
    print(f"g = {genus_sigma} + {genus_cap} = {g}")


solve_genus_problem()