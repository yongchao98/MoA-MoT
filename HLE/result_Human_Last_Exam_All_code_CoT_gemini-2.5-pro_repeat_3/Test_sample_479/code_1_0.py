def solve_genus_problem():
    """
    This function explains the solution to the maximal genus problem and
    provides a numerical illustration for the case of a torus (genus 1).
    """

    print("--- Mathematical Derivation ---")
    print("The problem is solved using a theorem by A. Ros in differential geometry.")
    print("The theorem states that any compact, connected, embedded surface that bounds a region in R^3 and has a genus of 1 or more, must have a point where the mean curvature is zero.")
    print("The problem states that the mean curvature *never* vanishes.")
    print("By the theorem's contrapositive, this implies the genus must be less than 1.")
    print("Since genus must be a non-negative integer, the only possibility is 0.")
    print("A sphere (genus 0) is a valid example, as it bounds a ball and has constant non-zero mean curvature (H=1/R).")
    print("Therefore, the maximal possible genus is 0.")
    
    print("\n--- Numerical Illustration: Why a Torus (Genus 1) Fails ---")
    print("Let's consider a torus as the boundary of a solid 'donut' region.")
    print("We analyze the mean curvature H at its outermost and innermost points, using the outward-pointing normal convention.")

    # Define the parameters for an example torus
    R = 4.0  # Major radius
    r = 1.0  # Minor radius
    
    print(f"\nTorus parameters: Major Radius R = {R}, Minor Radius r = {r}")
    
    # --- Calculation for the outermost point ---
    print("\n1. At the outermost point:")
    print("   Both principal curvatures are positive relative to the outward normal.")
    # The principal curvatures are k1 = 1/r and k2 = 1/(R+r)
    kappa1_out = 1 / r
    kappa2_out = 1 / (R + r)
    H_out = 0.5 * (kappa1_out + kappa2_out)
    
    print("   The equation for mean curvature is H = 0.5 * (k1 + k2).")
    print(f"   k1 = 1/r = 1/{r} = {kappa1_out:.3f}")
    print(f"   k2 = 1/(R+r) = 1/({R}+{r}) = {kappa2_out:.3f}")
    print(f"   H_outer = 0.5 * ({kappa1_out:.3f} + {kappa2_out:.3f}) = {H_out:.3f}")

    # --- Calculation for the innermost point ---
    print("\n2. At the innermost point:")
    print("   Both principal curvatures are negative relative to the outward normal.")
    # The principal curvatures are k1 = -1/r and k2 = -1/(R-r)
    kappa1_in = -1 / r
    kappa2_in = -1 / (R - r)
    H_in = 0.5 * (kappa1_in + kappa2_in)
    
    print("   The equation for mean curvature is H = 0.5 * (k1 + k2).")
    print(f"   k1 = -1/r = -1/{r} = {kappa1_in:.3f}")
    print(f"   k2 = -1/(R-r) = -1/({R}-{r}) = {kappa2_in:.3f}")
    print(f"   H_inner = 0.5 * ({kappa1_in:.3f} + {kappa2_in:.3f}) = {H_in:.3f}")
    
    # --- Conclusion from the illustration ---
    print("\n--- Conclusion ---")
    print(f"The mean curvature is positive ({H_out:.3f}) on the outside and negative ({H_in:.3f}) on the inside.")
    print("Since H is a continuous function on a connected surface, by the Intermediate Value Theorem, it must pass through 0 somewhere.")
    print("This confirms that a torus bounding a region cannot satisfy the condition of non-vanishing mean curvature.")
    print("\nThe maximal genus is 0.")

# Execute the function to print the full explanation and calculation.
solve_genus_problem()