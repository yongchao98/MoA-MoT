def solve_genus_problem():
    """
    Solves the topology problem by calculating the minimum genus for the containing surface.
    """
    # The given parameters for the initial surface, Σ.
    genus_sigma = 10
    boundary_components_sigma = 1

    print("--- Problem Analysis ---")
    print(f"We start with a surface Σ of genus h_Σ = {genus_sigma} and with {boundary_components_sigma} boundary component.")
    print("We want to find the smallest genus 'g' for a closed surface Σ' that contains Σ.")
    print("\nTo make Σ into a closed surface, we must attach a 'capping' surface, S_cap, to its boundary.")
    print("The genus of the final surface, g, is the sum of the genera of the two component surfaces.")
    print("The governing formula is: g = h_Σ + h_cap\n")

    print("--- Determining Minimum Genus ---")
    print("To find the minimum possible g, we must find the minimum possible genus for the capping surface, h_cap.")
    print("The problem states that the boundary of Σ is an 'unknotted' circle.")
    print("An unknotted circle in R^3 can be spanned by a smoothly embedded disk.")
    
    # A disk is a surface of genus 0 with one boundary component.
    min_genus_cap = 0
    
    print(f"A disk has a genus of 0. Therefore, the smallest possible genus for the capping surface is h_cap = {min_genus_cap}.")
    print("This choice of cap is always possible due to the 'unknotted' condition.\n")

    print("--- Final Calculation ---")
    print("We substitute the known values into our formula: g = h_Σ + h_cap")
    
    # Perform the final calculation
    min_genus_g = genus_sigma + min_genus_cap

    # Print the equation with the specific numbers, as requested.
    print(f"g = {genus_sigma} + {min_genus_cap}")
    print(f"g = {min_genus_g}\n")

    print(f"Thus, the smallest positive integer g such that Σ can be contained in a closed surface of genus g is {min_genus_g}.")

solve_genus_problem()