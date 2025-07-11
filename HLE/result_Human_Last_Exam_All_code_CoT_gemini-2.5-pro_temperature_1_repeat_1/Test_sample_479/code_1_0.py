def solve_maximal_genus():
    """
    This function explains the solution to find the maximal genus of a surface
    in R^3 with non-vanishing mean curvature.
    """
    
    print("Problem: What is the maximal genus of the boundary of a compact region in R^3, assuming its mean curvature never vanishes?")
    
    print("\n--- Step 1: Analysis of the condition ---")
    print("Let S be the surface. The condition is that its mean curvature H is never zero.")
    print("Since S is compact and connected, H must be either always positive (H > 0) or always negative (H < 0).")
    print("We can assume H > 0 everywhere without loss of generality.")
    
    print("\n--- Step 2: Checking possible genera ---")
    
    # Genus 0
    genus_0 = 0
    print(f"Case 1: Genus = {genus_0}")
    print("A sphere is a surface of genus 0. A sphere of radius R has constant mean curvature H = 1/R > 0.")
    print("Conclusion: Genus 0 is possible.")
    
    # Genus 1
    genus_1 = 1
    print(f"\nCase 2: Genus = {genus_1}")
    print("A torus is a surface of genus 1. A torus of revolution can have H > 0 everywhere if its major radius R is more than twice its minor radius r (R > 2r).")
    print("Conclusion: Genus 1 is possible.")

    # Genus >= 2
    print("\nCase 3: Genus >= 2")
    print("A fundamental theorem in differential geometry states that a compact, connected, embedded surface in R^3 with mean curvature of a fixed sign (like H > 0) can only be of genus 0 or 1.")
    print("This implies that surfaces of genus 2, 3, 4, etc., cannot have non-vanishing mean curvature.")
    
    print("\n--- Step 3: Final Conclusion ---")
    print("Based on the analysis, the possible genera are 0 and 1.")
    maximal_genus = 1
    print("The maximal genus is the largest possible value.")
    
    print("\nFinal Equation:")
    print(f"Maximal Genus = {maximal_genus}")

solve_maximal_genus()