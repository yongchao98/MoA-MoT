def solve_riemann_surface_problem():
    """
    Calculates the smallest and largest possible rank of the Néron-Severi group
    for the 15th symmetric power of a genus 3 Riemann surface.
    """

    # --- Problem Parameters ---
    g = 3  # Genus of the curve C
    d = 15 # The symmetric power

    # --- Introduction and Plan ---
    print("This script calculates the smallest and largest possible rank of the Néron-Severi group (NS) of X.")
    print(f"X is the {d}th symmetric power of a genus g = {g} Riemann surface C.")
    print("-" * 50)
    print("Plan:")
    print("1. Use the formula relating the rank of NS(X) to the rank of NS of the curve's Jacobian, J(C).")
    print("2. Determine the minimum and maximum possible ranks for NS(J(C)) for a genus 3 curve.")
    print("3. Compute the final minimum and maximum ranks for NS(X).")
    print("-" * 50)

    # --- Step 1: The Key Formula ---
    print("Step 1: The relationship between ρ(X) and ρ(J(C))")
    print("The rank of the Néron-Severi group of X = C^(d), denoted ρ(X), is related to the rank of the Néron-Severi group of the Jacobian J(C) by the formula:")
    print("ρ(C^(d)) = ρ(J(C)) + 1")
    print("This means we need to find the minimum and maximum of ρ(J(C)).")
    print("-" * 50)

    # --- Step 2: Bounds for ρ(J(C)) ---
    print("Step 2: Finding the bounds for ρ(J(C))")
    print(f"J(C) is an abelian variety of dimension g = {g}.")
    
    # Smallest rank for J(C)
    print("\nSmallest Rank (Generic Case):")
    print("For a generic curve, the Jacobian has the minimum possible rank for its Néron-Severi group.")
    min_rho_J = 1
    print(f"The minimum possible value is ρ_min(J(C)) = {min_rho_J}.")
    
    # Largest rank for J(C)
    print("\nLargest Rank (Complex Multiplication Case):")
    print(f"The rank of the Néron-Severi group of a g-dimensional abelian variety is at most g^2.")
    max_rho_J_bound = g**2
    print(f"For J(C) with g = {g}, this upper bound is {g}^2 = {max_rho_J_bound}.")
    print("This maximum is achieved for special curves (e.g., the Fermat quartic) whose Jacobians have maximal complex multiplication.")
    max_rho_J = max_rho_J_bound
    print(f"The maximum possible value is ρ_max(J(C)) = {max_rho_J}.")
    print("-" * 50)

    # --- Step 3: Final Calculation ---
    print("Step 3: Calculating the smallest and largest ranks for X = C^(15)")
    
    # Smallest rank for X
    min_rho_X = min_rho_J + 1
    print(f"\nThe smallest possible rank of NS(X) is:")
    print(f"ρ_min(X) = ρ_min(J(C)) + 1 = {min_rho_J} + 1 = {min_rho_X}")

    # Largest rank for X
    max_rho_X = max_rho_J + 1
    print(f"\nThe largest possible rank of NS(X) is:")
    print(f"ρ_max(X) = ρ_max(J(C)) + 1 = {max_rho_J} + 1 = {max_rho_X}")

# Execute the function to print the solution
solve_riemann_surface_problem()