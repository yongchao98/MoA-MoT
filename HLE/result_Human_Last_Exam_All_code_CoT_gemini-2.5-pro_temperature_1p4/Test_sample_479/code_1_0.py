import math

def find_maximal_genus():
    """
    This function determines the maximal genus of a compact region's boundary
    in R^3 with non-vanishing mean curvature, based on the Gauss-Bonnet theorem
    and a related geometric inequality.
    """
    print("To find the maximal genus, we follow these steps:")
    print("1. Recall the Gauss-Bonnet Theorem for a compact surface S with genus g:")
    print("   Integral(K dA) = 4 * pi * (1 - g)")
    print("   where K is the Gauss curvature.")
    
    print("\n2. Use the condition that the mean curvature H is non-vanishing (e.g., H > 0 everywhere).")
    print("   A known geometric inequality states that for such a surface, the total Gauss curvature must be non-negative:")
    print("   Integral(K dA) >= 0")

    print("\n3. Combine the theorem and the inequality:")
    print("   4 * pi * (1 - g) >= 0")
    print("   This simplifies to the condition: g <= 1")
    
    print("\nThis inequality implies the genus can be at most 1. Let's verify for integer genera:")
    
    max_g_found = -1
    
    # We check genera 0, 1, and 2 to illustrate the principle.
    for g in range(3):
        # According to Gauss-Bonnet, the integral of K depends on g.
        # Equation: integral_K = 4 * pi * (1 - g)
        integral_K = 4 * math.pi * (1 - g)
        
        # Check if this value is consistent with the inequality Integral(K dA) >= 0
        is_possible = (integral_K >= 0)
        
        print(f"\n- For genus g = {g}:")
        
        # Print the full equation for this genus
        one_minus_g = 1 - g
        print(f"  The equation is Integral(K dA) = 4 * pi * (1 - {g}) = 4 * pi * ({one_minus_g}) = {integral_K:.4f}")

        if is_possible:
            print(f"  This is POSSIBLE because the result ({integral_K:.4f}) is non-negative.")
            if g > max_g_found:
                max_g_found = g
        else:
            print(f"  This is IMPOSSIBLE because the result ({integral_K:.4f}) is negative, violating the H>0 condition.")
    
    print("\n4. Confirm that genus 1 is achievable.")
    print("   A torus (genus g=1) can be constructed with mean curvature H > 0 everywhere.")
    
    print("\nConclusion: The analysis shows that the genus cannot exceed 1, and g=1 is achievable.")
    print(f"Therefore, the maximal genus is {max_g_found}.")

# Execute the function to print the analysis.
find_maximal_genus()