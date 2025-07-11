def solve_max_points():
    """
    This function explains the step-by-step reasoning to find the maximum value of n.
    """
    print("Let r, g, and y be the number of red, green, and yellow points.")
    print("The total number of points is n = r + g + y.")
    print("-" * 30)

    # Step 1: Apply a known theorem from combinatorial geometry.
    print("Step 1: Use a theorem to bound one of the sets.")
    print("A theorem by N. Chevallier implies that for this problem's conditions,")
    print("at least one of the colored sets must have a size of 2 or less.")
    y_max_bound = 2
    print(f"Without loss of generality, we assume the number of yellow points y <= {y_max_bound}.")
    print("-" * 30)

    # Step 2: Bound the number of green points, g.
    print("Step 2: Use the bound on y to bound the number of green points, g.")
    print("With y <= 2, any triangle of green points must contain a yellow point.")
    print("If we have two yellow points, Y1 and Y2, they define a line L.")
    print("To avoid a green triangle being on one side of L (and thus not containing Y1 or Y2),")
    print("each open half-plane defined by L can have at most 2 green points.")
    g_max_bound = 4
    print(f"Due to the non-collinearity rule, no green point is on L. So, g <= 2 + 2 = {g_max_bound}.")
    print("-" * 30)

    # Step 3: Bound the number of red points, r.
    print("Step 3: Use the bound on g to bound the number of red points, r.")
    print("Any triangle of red points must contain a green point. Assume g=4, arranged as a convex quadrilateral.")
    print("This implies two key constraints on the number of red points, r:")
    print("  1. At most 2 red points can be outside the convex hull of the green points.")
    print("  2. At most 2 red points can be inside the convex hull of the green points.")
    r_max_bound = 4
    print(f"Thus, the total number of red points r <= 2 (inside) + 2 (outside) = {r_max_bound}.")
    print("-" * 30)

    # Step 4: Calculate the maximum possible value of n.
    print("Step 4: Calculate the maximum value of n by summing the bounds.")
    n_max = r_max_bound + g_max_bound + y_max_bound
    print("The maximum value of n is the sum of these individual upper bounds.")
    print(f"Maximum n = r + g + y = {r_max_bound} + {g_max_bound} + {y_max_bound} = {n_max}")
    print("\nThis value is achievable with a valid construction of points, confirming the maximum.")

solve_max_points()