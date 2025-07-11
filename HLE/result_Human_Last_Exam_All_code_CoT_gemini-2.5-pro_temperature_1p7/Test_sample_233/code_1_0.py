def solve_surface_genus():
    """
    Calculates the smallest positive integer g for the problem statement.
    """
    # The genus of the initial smoothly embedded oriented surface Sigma.
    genus_sigma = 10
    
    print("Step 1: Understand the Goal")
    print("We start with a surface Sigma of genus 10 with one boundary.")
    print("We need to find the smallest genus 'g' of a closed surface Sigma' that contains Sigma, valid for any such Sigma.")
    print("-" * 20)

    print("Step 2: Forming the Closed Surface")
    print("A closed surface has no boundary. To create Sigma', we must 'cap' the boundary of Sigma.")
    print("This is done by gluing a capping surface, let's call it F, to the boundary of Sigma.")
    print("-" * 20)

    print("Step 3: The Genus Formula")
    print("The genus of the resulting surface Sigma' is the sum of the genera of the original surface Sigma and the capping surface F.")
    print("The formula is: g(Sigma') = g(Sigma) + g(F)")
    print("-" * 20)
    
    print("Step 4: The Worst-Case Scenario")
    print("The problem requires a solution for 'regardless of our choice of Sigma'. This means we must consider the worst-case embedding of Sigma in 3D space.")
    print("A 'handle' of the surface Sigma can be topologically linked with its own boundary. To create an embedded closing surface Sigma', the cap F must not intersect Sigma.")
    print("To avoid intersection with a linked handle, the cap F must have its own handle.")
    print(f"The number of handles on Sigma is equal to its genus, which is {genus_sigma}.")
    print(f"In the worst-case configuration, all {genus_sigma} handles of Sigma are linked with its boundary.")
    
    # The genus of the capping surface F in the worst-case scenario
    # must be equal to the number of linked handles it needs to avoid.
    genus_f_worst_case = genus_sigma
    print(f"Therefore, the capping surface F must have a genus of {genus_f_worst_case} to avoid all the linked handles.")
    print("-" * 20)

    print("Step 5: Final Calculation")
    # The final genus g is the genus of Sigma' in the worst-case scenario.
    final_genus = genus_sigma + genus_f_worst_case
    
    print(f"We calculate the smallest possible genus g for Sigma' that works in all cases:")
    print(f"g = g(Sigma) + g(F_worst_case)")
    print(f"g = {genus_sigma} + {genus_f_worst_case} = {final_genus}")
    
    print("-" * 20)
    print(f"The smallest positive integer g is {final_genus}.")

solve_surface_genus()
<<<20>>>