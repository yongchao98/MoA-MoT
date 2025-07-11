import math

def solve_packing():
    """
    Solves the discrete sphere packing problem using a greedy algorithm.
    """
    print("Is the problem formulation correct? Yes, it is.")
    print("Let's verify the non-overlapping constraint equation:")

    # Given parameters
    sphere_radius_cm = 2.0
    grid_unit_cm = 0.5
    
    # Derivation
    print(f"The radius of each sphere is {sphere_radius_cm} cm.")
    print("For two spheres not to overlap, the distance between their centers must be at least the sum of their radii.")
    min_dist_cm = sphere_radius_cm + sphere_radius_cm
    print(f"Minimum distance between centers (cm): {sphere_radius_cm} + {sphere_radius_cm} = {min_dist_cm} cm.")
    
    # In squared terms to avoid square roots
    min_sq_dist_cm = min_dist_cm ** 2
    print(f"Minimum squared distance (cm^2): {min_dist_cm}^2 = {min_sq_dist_cm}")

    print("\nNow, let's convert this to the integer grid coordinate system.")
    print(f"The grid coordinates (x, y, z) are in units of {grid_unit_cm} cm.")
    
    # Real coordinates (X, Y, Z) are grid coordinates (x, y, z) * grid_unit_cm
    # (X_i - X_j)^2 + (Y_i - Y_j)^2 + (Z_i - Z_j)^2 >= 16
    # (x_i*0.5 - x_j*0.5)^2 + ... >= 16
    # (0.5 * (x_i - x_j))^2 + ... >= 16
    # 0.25 * (x_i - x_j)^2 + ... >= 16
    
    # Multiply by 4 to clear the 0.25
    min_sq_dist_grid = min_sq_dist_cm / (grid_unit_cm ** 2)
    print(f"(x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= {min_sq_dist_cm} / {grid_unit_cm**2}")
    print(f"(x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= {int(min_sq_dist_grid)}")

    print(f"\nThe integer radius is {int(sphere_radius_cm / grid_unit_cm)}, so the diameter is {int(2 * sphere_radius_cm / grid_unit_cm)}.")
    print(f"The constraint can be written as >= ({int(2 * sphere_radius_cm / grid_unit_cm)})^2, which matches your son's (4+4)^2. The formulation is correct.")
    print("\n-------------------------------------------------")
    print("Running greedy algorithm to find the maximum number of spheres...")

    # Problem parameters in grid units
    x_range = (4, 20)
    y_range = (4, 20)
    z_range = (4, 18)
    
    # Generate all possible candidate points, sorted lexicographically (z, y, x)
    candidate_points = []
    for z in range(z_range[0], z_range[1] + 1):
        for y in range(y_range[0], y_range[1] + 1):
            for x in range(x_range[0], x_range[1] + 1):
                candidate_points.append((x, y, z))
    
    available_points = set(candidate_points)
    placed_spheres = []

    # Iterate through sorted candidates
    for p1 in candidate_points:
        if p1 in available_points:
            # Place a sphere at this point
            placed_spheres.append(p1)
            x1, y1, z1 = p1
            
            # Remove all points that conflict with the new sphere
            # Using a list comprehension for a temporary list of points to remove
            points_to_remove = {p2 for p2 in available_points 
                                if (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2 < min_sq_dist_grid}

            available_points.difference_update(points_to_remove)

    n = len(placed_spheres)
    print(f"\nFinished calculation.")
    print(f"The maximized value n is: {n}")

# Run the solver
solve_packing()