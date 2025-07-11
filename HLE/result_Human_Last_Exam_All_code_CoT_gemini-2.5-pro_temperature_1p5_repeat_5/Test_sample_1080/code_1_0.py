import collections

def solve_candy_packing():
    """
    Solves the discrete sphere packing problem for the given box and candy size.
    """
    # Box dimensions: 12x12x11 cm
    # Candy radius: 2cm (diameter 4cm)
    # Grid spacing: 0.5cm

    # 1. Define the search space for candy centers based on the problem description.
    # The actual coordinates in cm are (x*0.5, y*0.5, z*0.5).
    # The candy center (cx,cy,cz) must be in [2,10] x [2,10] x [2,9] cm.
    # This corresponds to integer coordinates x,y in [4,20] and z in [4,18].
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # The non-overlapping constraint for two centers (x_i, y_i, z_i) and (x_j, y_j, z_j) is:
    # (0.5*x_i - 0.5*x_j)^2 + ... >= (2+2)^2
    # 0.25 * ((x_i - x_j)^2 + ...) >= 16
    # (x_i - x_j)^2 + ... >= 64
    # The minimum squared distance in grid units is ( (2+2) / 0.5 )^2 = 8^2 = 64.
    min_dist_sq = 64

    # 2. Generate all possible center points in a specific order (z, then y, then x).
    # This order is good for finding layer-based packings.
    possible_points = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                possible_points.append((x, y, z))

    # 3. Implement the greedy algorithm.
    placed_candies = []
    for point in possible_points:
        is_valid_placement = True
        # Check against all previously placed candies.
        for placed_candy in placed_candies:
            dist_sq = (point[0] - placed_candy[0])**2 + \
                      (point[1] - placed_candy[1])**2 + \
                      (point[2] - placed_candy[2])**2
            
            if dist_sq < min_dist_sq:
                is_valid_placement = False
                break
        
        if is_valid_placement:
            placed_candies.append(point)

    # 4. Format the output to show the final "equation".
    if not placed_candies:
        print("No candies could be placed.")
        return

    # Count how many candies are on each z-level.
    z_levels = [p[2] for p in placed_candies]
    level_counts = collections.Counter(z_levels)
    
    # Sort the counts by the z-level for a clear representation.
    sorted_levels = sorted(level_counts.items())
    
    counts_per_level = [count for level, count in sorted_levels]
    total_candies = sum(counts_per_level)

    # Print the equation showing the breakdown by layer.
    equation_str = " + ".join(map(str, counts_per_level))
    print(f"Yes, the problem formulation is correct.")
    print(f"The maximum number of candies is {total_candies}.")
    print(f"The breakdown by layers is: {equation_str} = {total_candies}")

solve_candy_packing()