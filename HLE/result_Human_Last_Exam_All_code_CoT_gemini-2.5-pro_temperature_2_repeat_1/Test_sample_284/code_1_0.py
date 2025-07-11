def solve_fourier_restriction_problem():
    """
    Calculates the critical exponent for the Fourier restriction problem on the moment curve in R^3.
    """
    # The vertices of the sharp boundedness region for the extension operator in the (1/p, 1/q) plane
    # are the origin and the points P1, P2, P3.
    p1 = (1/6, 1/2)
    p2 = (1/8, 1/4)
    p3 = (1/12, 0)
    
    # We want to find the range of p for which an L^p function with Fourier support on the curve can exist.
    # This corresponds to the projection of the boundedness region onto the 1/p axis.
    one_over_p_coords = [p1[0], p2[0], p3[0], 0]
    
    # The minimal value of p corresponds to the maximal value of 1/p.
    max_one_over_p = max(one_over_p_coords)
    
    # The critical exponent p_c is the inverse of this maximal value.
    p_critical = 1 / max_one_over_p
    
    print("The vertices defining the boundedness region in the (1/p, 1/q) plane are:")
    print(f"(0, 0), {p1}, {p2}, {p3}")
    
    print("\nThe coordinates for 1/p at these vertices are:")
    # Print each number used in the calculation.
    print(f"0, {p1[0]:.4f} (from 1/6), {p2[0]:.4f} (from 1/8), {p3[0]:.4f} (from 1/12)")
    
    print(f"\nThe maximum value of 1/p in this set is 1/6.")
    
    # Show the final calculation
    print(f"\nThe largest possible value of p for which NO L^p function can have its Fourier support on the moment curve is the critical value p_c.")
    print(f"p_c = 1 / max(0, 1/6, 1/8, 1/12) = 1 / (1/6) = {int(p_critical)}")

solve_fourier_restriction_problem()
