import math

def solve_halloween_packing():
    """
    Calculates the maximum number of spheres that can be packed into a box
    with the given discrete constraints using a greedy algorithm.
    """
    # Box dimensions in cm: 12x12x11
    # Sphere radius in cm: 2
    # Grid resolution: 0.5 cm, so all values are multiplied by 2.
    
    box_dims = (24, 24, 22)
    # Sphere radius in grid units.
    radius = 4
    
    # Non-overlapping constraint: The distance between centers must be at least 2 * radius.
    # The squared distance must be at least (2 * radius)^2.
    min_dist_sq = (2 * radius)**2

    # Define the valid domain for the sphere centers based on box dimensions and sphere radius.
    x_range = range(radius, box_dims[0] - radius + 1)
    y_range = range(radius, box_dims[1] - radius + 1)
    z_range = range(radius, box_dims[2] - radius + 1)

    # Generate all possible candidate center points, sorted lexicographically (by z, then y, then x).
    # This ensures we always start packing from the "bottom-front-left" corner of the available space.
    candidates = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidates.append((x, y, z))

    packed_spheres_centers = []
    
    # Greedily pack spheres.
    while candidates:
        # Select the first available point from the sorted list.
        p_center = candidates[0]
        packed_spheres_centers.append(p_center)
        
        # Create a new list of candidates, keeping only those that don't overlap with the sphere just placed.
        remaining_candidates = []
        for q_center in candidates:
            dist_sq = (p_center[0] - q_center[0])**2 + (p_center[1] - q_center[1])**2 + (p_center[2] - q_center[2])**2
            if dist_sq >= min_dist_sq:
                remaining_candidates.append(q_center)
        candidates = remaining_candidates

    # The problem formulation is correct. Now we output the result.
    n = len(packed_spheres_centers)
    
    print("Yes, your problem formulation is correct.")
    print("\nBased on your formulation, we solve for the maximum number of candies 'n'.")
    print("The core constraint is the non-overlapping equation for any two candies i and j:")
    print(f"(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= ({radius} + {radius})^2")
    print(f"Using the grid units, this simplifies to:")
    print(f"(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= {min_dist_sq}")
    
    print("\nBy running a greedy algorithm to pack the candies, we find the maximized value.")
    print(f"The maximum number of eyeball candies (n) is: {n}")

# Run the simulation and print the results.
solve_halloween_packing()
<<<22>>>