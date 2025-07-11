def solve_genus_problem():
    """
    Calculates the smallest possible genus g for a closed surface containing a given surface Sigma.
    """
    # Step 1: Define properties of the initial surface Sigma
    genus_sigma = 10
    boundaries_sigma = 1
    print(f"The initial surface Sigma has genus h = {genus_sigma} and b = {boundaries_sigma} boundary component.")

    # Calculate the Euler characteristic of Sigma
    chi_sigma = 2 - 2 * genus_sigma - boundaries_sigma
    print(f"The Euler characteristic of Sigma is calculated as chi(Sigma) = 2 - 2*h - b.")
    print(f"chi(Sigma) = 2 - 2*({genus_sigma}) - {boundaries_sigma} = {chi_sigma}")
    print("-" * 20)

    # Step 2: Define the cap to close the boundary
    # To find the smallest resulting genus g, we must use the simplest possible cap.
    # The simplest surface with one boundary component is a disk, which has genus 0.
    genus_cap = 0
    boundaries_cap = 1
    print(f"To close the boundary, we attach a cap. To get the smallest final genus, we use the simplest cap: a disk.")
    print(f"The disk has genus h_cap = {genus_cap} and b_cap = {boundaries_cap} boundary.")
    
    # Calculate the Euler characteristic of the disk (cap)
    chi_cap = 2 - 2 * genus_cap - boundaries_cap
    print(f"The Euler characteristic of the disk is chi(Cap) = 2 - 2*({genus_cap}) - {boundaries_cap} = {chi_cap}")
    print("-" * 20)

    # Step 3: Calculate the Euler characteristic of the final closed surface Sigma'
    # The intersection is the boundary circle, which has an Euler characteristic of 0.
    chi_intersection = 0
    print(f"The new surface Sigma' is formed by gluing Sigma and the disk cap.")
    print(f"Its Euler characteristic chi(Sigma') is chi(Sigma) + chi(Cap) - chi(Boundary).")
    
    chi_sigma_prime = chi_sigma + chi_cap - chi_intersection
    print(f"chi(Sigma') = {chi_sigma} + {chi_cap} - {chi_intersection} = {chi_sigma_prime}")
    print("-" * 20)
    
    # Step 4: Calculate the genus g of the new closed surface Sigma'
    # For a closed surface, chi = 2 - 2g. So, g = (2 - chi) / 2.
    print(f"For a closed surface like Sigma', the genus g is related to its Euler characteristic by chi = 2 - 2g.")
    print(f"Solving for g: g = (2 - chi(Sigma')) / 2")
    
    # We use integer division // as genus must be an integer
    g = (2 - chi_sigma_prime) // 2
    
    print(f"g = (2 - ({chi_sigma_prime})) / 2 = {2 - chi_sigma_prime} / 2 = {g}")
    print("-" * 20)
    
    print(f"The smallest positive integer g is {g}.")

solve_genus_problem()
