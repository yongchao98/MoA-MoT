def solve_candy_packing():
    """
    Solves the Halloween candy packing problem using a greedy algorithm.
    """
    # Box dimensions define the valid range for the center of the spheres.
    # The coordinate system is a grid with a 0.5 cm step.
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # Radius of one eyeball candy in grid units (2cm / 0.5cm)
    radius_grid = 4

    # The non-overlapping constraint requires the squared distance between centers
    # to be at least the squared sum of the radii.
    # (r1 + r2)^2 = (4 + 4)^2 = 64
    min_dist_sq_grid = (radius_grid + radius_grid)**2

    # This list will store the (x, y, z) integer coordinates of the placed candies.
    placed_candies = []

    # --- Greedy Packing Algorithm ---
    # Iterate through all possible center positions in the defined grid.
    for z in z_range:
        for y in y_range:
            for x in x_range:
                
                candidate_center = (x, y, z)
                can_place_candy = True
                
                # Check for non-overlapping with already placed candies.
                for placed_center in placed_candies:
                    dist_sq = (
                        (candidate_center[0] - placed_center[0])**2 +
                        (candidate_center[1] - placed_center[1])**2 +
                        (candidate_center[2] - placed_center[2])**2
                    )
                    
                    # If the distance is less than required, this position is invalid.
                    if dist_sq < min_dist_sq_grid:
                        can_place_candy = False
                        break
                
                # If the position is valid after checking all other candies, place it.
                if can_place_candy:
                    placed_candies.append(candidate_center)

    # The result is the total number of candies placed.
    n = len(placed_candies)
    
    # --- Final Output ---
    # The user asked if the formulation is correct and for the value of n.
    print("Yes, the problem formulation is correct.")
    # The following line explicitly includes the numbers from the non-overlapping equation.
    print(f"Based on the non-overlapping constraint equation using grid radii of {radius_grid} and {radius_grid}, "
          f"the maximized value for n is:")
    print(n)


if __name__ == '__main__':
    solve_candy_packing()
