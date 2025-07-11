def solve_genus_problem():
    """
    This function explains the reasoning to find the maximal genus of the described surface.
    """
    
    print("Step 1: The problem describes a compact region Omega whose boundary S is a smooth, connected surface.")
    print("The key condition is that its mean curvature H is never zero.")
    
    print("\nStep 2: By the maximum principle, it can be shown that if S is the boundary of a compact region, its mean curvature (with respect to the outward normal) must be negative everywhere, H < 0.")
    
    print("\nStep 3: If we consider the inward normal, the mean curvature becomes positive everywhere, H > 0.")
    
    print("\nStep 4: The Heinz-Hopf theorem states that for a compact surface in R^3 with H > 0, the degree of its Gauss map, deg(G), must be non-negative.")
    print("So, we have the condition: deg(G) >= 0.")
    
    print("\nStep 5: The Gauss-Bonnet theorem relates the degree of the Gauss map to the genus (g) of the surface:")
    print("deg(G) = 1 - g.")
    
    print("\nStep 6: Combining these results gives the final inequality.")
    
    # The numbers in the final equation
    one = 1
    zero = 0
    
    print(f"The inequality is: deg(G) = {one} - g >= {zero}.")
    print("Solving for g, we get: g <= 1.")
    
    print("\nStep 7: We know that genus 0 (a sphere) and genus 1 (a suitable torus) are possible examples.")
    
    max_genus = 1
    print(f"\nConclusion: The maximal possible genus for the boundary of such a region is {max_genus}.")

solve_genus_problem()