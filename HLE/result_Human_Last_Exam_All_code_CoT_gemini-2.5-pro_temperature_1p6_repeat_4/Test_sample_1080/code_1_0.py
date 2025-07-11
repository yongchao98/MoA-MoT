import collections

def solve_candy_packing():
    """
    Solves the sphere packing problem using a greedy algorithm.
    """
    # Define the boundaries for the center coordinates in the 0.5cm grid
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # The non-overlapping constraint: distance between centers must be >= 8 units.
    # We use the squared distance to avoid square roots.
    min_dist_sq = (4 + 4)**2  # (radius_grid + radius_grid)^2 = 8^2 = 64

    # Generate a list of all possible candidate center points.
    # The order (z, then y, then x) makes the greedy algorithm fill bottom-up.
    candidate_points = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_points.append((x, y, z))

    # This list will store the coordinates of the candies we pack.
    packed_spheres = []

    # Iterate through each candidate point and try to place a candy.
    for p1 in candidate_points:
        is_valid_placement = True
        # Check against all candies that have already been placed.
        for p2 in packed_spheres:
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            if dist_sq < min_dist_sq:
                is_valid_placement = False
                break
        
        if is_valid_placement:
            packed_spheres.append(p1)

    # Now we analyze the packing to present the result.
    print("The problem formulation is correct.\n")
    print("Based on the greedy packing algorithm, here is the result:\n")

    if not packed_spheres:
        print("No candies could be packed.")
        n = 0
    else:
        # Count how many spheres were placed in each z-layer.
        counts_per_z_layer = collections.Counter(z for x, y, z in packed_spheres)
        sorted_layers = sorted(counts_per_z_layer.items())

        sum_parts = []
        # Print the breakdown per layer.
        for z_coord, count in sorted_layers:
            # Convert grid z-coordinate back to cm for clarity.
            real_z = z_coord * 0.5
            print(f"Layer with centers at z = {real_z:.1f} cm (grid z={z_coord}): {count} candies")
            sum_parts.append(str(count))
        
        n = len(packed_spheres)
        equation_str = " + ".join(sum_parts)
        
        print("\nThe final calculation for the total number of candies is:")
        print(f"Total n = {equation_str} = {n}")

solve_candy_packing()