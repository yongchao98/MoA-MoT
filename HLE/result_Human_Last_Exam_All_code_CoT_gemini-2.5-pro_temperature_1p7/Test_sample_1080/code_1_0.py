import math

def solve_halloween_packing():
    """
    Solves the integer-grid sphere packing problem using a greedy algorithm.
    """
    # Define the valid range for the integer coordinates of the sphere centers.
    x_min, x_max = 4, 20
    y_min, y_max = 4, 20
    z_min, z_max = 4, 18

    # The squared distance between centers must be at least (2*radius/grid_step)^2
    # (2 * 2 cm / 0.5 cm)^2 = (4 / 0.5)^2 = 8^2 = 64.
    min_sq_dist = 64

    # Generate all possible center locations and keep them as a set for efficient removal.
    available_points = set()
    for z in range(z_min, z_max + 1):
        for y in range(y_min, y_max + 1):
            for x in range(x_min, x_max + 1):
                available_points.add((x, y, z))

    # Create a deterministically sorted list of candidate points to iterate through.
    # Sorting by z, then y, then x ensures we fill layer by layer from a corner.
    candidate_points = sorted(list(available_points), key=lambda p: (p[2], p[1], p[0]))
    
    selected_centers = []

    # Greedily select points.
    for candidate in candidate_points:
        # Check if the candidate point is still available.
        if candidate in available_points:
            # If so, select this point as a new center.
            selected_centers.append(candidate)
            
            # Now, remove all points from the available set that conflict with this new center.
            # This includes the candidate point itself.
            cx, cy, cz = candidate
            to_remove = set()

            # For efficiency, only check points in a bounding box around the candidate.
            # The radius of conflict is sqrt(64) = 8.
            z_check_range = range(max(z_min, cz - 7), min(z_max + 1, cz + 8))
            y_check_range = range(max(y_min, cy - 7), min(y_max + 1, cy + 8))
            x_check_range = range(max(x_min, cx - 7), min(x_max + 1, cx + 8))

            for z in z_check_range:
                for y in y_check_range:
                    for x in x_check_range:
                        p2 = (x, y, z)
                        # We only need to check points that are actually still available.
                        if p2 in available_points:
                            dist_sq = (cx - x)**2 + (cy - y)**2 + (cz - z)**2
                            if dist_sq < min_sq_dist:
                                to_remove.add(p2)
            
            available_points -= to_remove

    # The "final equation" is simply the number of spheres found.
    print(f"The maximized number of eyeball candies is: {len(selected_centers)}")
    
    # Optional: print the coordinates of each candy
    # print("The coordinates of the centers are:")
    # for center in sorted(selected_centers):
    #     print(center)


if __name__ == '__main__':
    solve_halloween_packing()
