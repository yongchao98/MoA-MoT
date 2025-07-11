def solve_genus_problem():
    """
    This function explains the logical steps to find the maximal genus
    of a surface with non-vanishing mean curvature.
    """

    # The reasoning leads to an inequality for the genus g.
    # From a theorem in differential geometry, if Mean Curvature H >= 0, then the
    # integral of Gaussian Curvature K dA >= 0.
    # The Gauss-Bonnet theorem states: integral(K dA) = 4 * pi * (1 - g).
    # Combining these, we get the inequality for g.

    inequality_1 = "4 * pi * (1 - g) >= 0"
    inequality_2 = "1 - g >= 0"
    inequality_3 = "g <= 1"
    
    # We know that genus g=0 (sphere) and g=1 (torus) surfaces with H>0 exist.
    max_genus = 1

    print("Step 1: From the properties of the surface and its mean curvature, we deduce that the mean curvature H must be strictly positive everywhere (H > 0), by choosing the correct orientation.")
    print("Step 2: A theorem in differential geometry states that for such a surface, the total Gaussian curvature must be non-negative.")
    print("Step 3: The Gauss-Bonnet Theorem relates the total Gaussian curvature to the genus g.")
    print("Step 4: Combining these results gives the inequality:")
    print(f"    {inequality_1}")
    print("Step 5: Since 4*pi is positive, the inequality simplifies to:")
    print(f"    {inequality_2}")
    print("Step 6: This gives a final constraint on the genus g:")
    print(f"    {inequality_3}")
    print(f"Step 7: We know that a surface of genus {max_genus} with non-vanishing mean curvature exists (a specific type of torus).")
    print(f"Therefore, the maximal possible genus is {max_genus}.")

solve_genus_problem()