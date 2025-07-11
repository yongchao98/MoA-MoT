import collections

def solve_packing():
    """
    Solves the sphere packing problem using a greedy algorithm.
    """
    # Box dimensions are 12x12x11 cm. Sphere radius is 2cm.
    # Grid step is 0.5 cm.
    # Sphere radius in grid units is 2 / 0.5 = 4.
    # Diameter in grid units is 8.
    # The minimum squared distance between centers is (diameter)^2 = 8*8 = 64.
    min_sq_dist = 64

    # The valid integer ranges for sphere centers (x, y, z).
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # List to store the coordinates of the spheres we place.
    placed_spheres = []

    # Iterate through all possible points in a specific order (z, then y, then x).
    # This corresponds to filling the box layer by layer from the bottom up.
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_point = (x, y, z)
                is_valid = True
                # Check for collision with already placed spheres.
                for placed_point in placed_spheres:
                    # Calculate squared distance. It's faster than sqrt.
                    dist_sq = (candidate_point[0] - placed_point[0])**2 + \
                              (candidate_point[1] - placed_point[1])**2 + \
                              (candidate_point[2] - placed_point[2])**2
                    
                    if dist_sq < min_sq_dist:
                        is_valid = False
                        break # Collision found, no need to check further.
                
                if is_valid:
                    placed_spheres.append(candidate_point)

    # The maximum number of spheres is the total count.
    n = len(placed_spheres)

    # To create the equation, we count how many spheres are in each z-layer.
    layer_counts = collections.defaultdict(int)
    for sphere in placed_spheres:
        layer_counts[sphere[2]] += 1

    # Sort the layers by z-coordinate to get a sensible order for the equation.
    sorted_layers = sorted(layer_counts.items())
    
    # Format the numbers for the final equation string.
    count_strings = [str(count) for z, count in sorted_layers]
    equation_str = " + ".join(count_strings)

    print(f"Yes, the problem formulation is correct.")
    print(f"The maximum number of eyeball candies is: {equation_str} = {n}")

solve_packing()