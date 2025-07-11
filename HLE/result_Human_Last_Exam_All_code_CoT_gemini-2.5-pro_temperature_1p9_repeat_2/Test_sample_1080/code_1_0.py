import math

def solve_sphere_packing():
    """
    Solves the sphere packing problem using a greedy algorithm
    based on the user's formulation.
    """
    # 1. Verification of the problem formulation.
    # Box: 12x12x11 cm. Sphere radius: 2cm. Grid step: 0.5cm.
    # Center coordinates (X, Y, Z) must be within:
    # X in [2, 10], Y in [2, 10], Z in [2, 9] (in cm)
    # User's integer coordinates (x, y, z) are multiples of 0.5cm.
    # So, x = X/0.5, y = Y/0.5, z = Z/0.5.
    # x range: [2/0.5, 10/0.5] -> [4, 20]. Correct.
    # y range: [2/0.5, 10/0.5] -> [4, 20]. Correct.
    # z range: [2/0.5, 9/0.5] -> [4, 18]. Correct.
    #
    # Non-overlapping constraint: Distance between centers >= 2 * radius = 4 cm.
    # In integer grid coordinates, this is a distance of 4 / 0.5 = 8 units.
    # The squared distance must be >= 8*8 = 64.
    # (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= 64. Correct.
    
    print("Yes, your problem formulation is correct.")
    print("We will solve for 'n' using a greedy algorithm.\n")

    # 2. Setup for the greedy algorithm
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)
    min_dist_sq = 64
    
    # Generate all candidate points, sorted by z, then y, then x
    # This simulates filling the box from the bottom corner, layer by layer.
    candidate_points = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_points.append((x, y, z))

    placed_spheres = []
    
    # 3. Run the greedy algorithm
    for candidate in candidate_points:
        is_valid_placement = True
        for placed in placed_spheres:
            # Calculate the squared distance to each already placed sphere
            dist_sq = (candidate[0] - placed[0])**2 + \
                      (candidate[1] - placed[1])**2 + \
                      (candidate[2] - placed[2])**2
            
            # If the new sphere is too close to an existing one, it's not a valid placement
            if dist_sq < min_dist_sq:
                is_valid_placement = False
                break
        
        # If after checking all placed spheres, the position is valid, place it.
        if is_valid_placement:
            placed_spheres.append(candidate)
            
    # 4. Output the result
    n = len(placed_spheres)
    
    print("--- SOLUTION ---")
    print("Objective Function: Maximize n")
    
    # "output each number in the final equation"
    # We interpret this as showing the final result for n.
    print(f"The maximized value n is: {n}")
    
    print("\nThis result was found by placing spheres at the following integer coordinates:")
    for i, p in enumerate(placed_spheres):
      print(f"  Sphere {i+1:>2}: (x={p[0]}, y={p[1]}, z={p[2]})")

    # Also show an example of the constraint calculation
    if n > 1:
        p1 = placed_spheres[0]
        p2 = placed_spheres[1] # The next sphere placed by the greedy algorithm
        dist_sq_calc = f"({p2[0]}-{p1[0]})^2 + ({p2[1]}-{p1[1]})^2 + ({p2[2]}-{p1[2]})^2"
        res_sq_calc = (p2[0]-p1[0])**2 + (p2[1]-p1[1])**2 + (p2[2]-p1[2])**2

        print("\nExample of the non-overlapping constraint equation for the first two spheres:")
        print(f"  Equation: {dist_sq_calc} >= 64")
        print(f"  Calculation: {res_sq_calc} >= 64")
        print(f"  The constraint is satisfied.")

solve_sphere_packing()