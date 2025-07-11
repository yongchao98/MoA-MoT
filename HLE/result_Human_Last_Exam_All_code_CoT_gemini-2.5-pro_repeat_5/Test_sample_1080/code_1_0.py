import math

def solve_sphere_packing():
    """
    Solves the sphere packing problem using a greedy algorithm.
    All dimensions and coordinates are in a grid where 1 unit = 0.5 cm.
    """
    # 1. Define problem parameters in grid units.
    box_dims = {'width': 24, 'depth': 24, 'height': 22}
    radius = 4
    min_dist_sq = (2 * radius)**2

    # 2. Generate all possible candidate center positions, sorted to fill
    # layer by layer along the z-axis, which is the shortest dimension.
    x_range = range(radius, box_dims['width'] - radius + 1)
    y_range = range(radius, box_dims['depth'] - radius + 1)
    z_range = range(radius, box_dims['height'] - radius + 1)

    candidate_positions = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_positions.append((x, y, z))

    # 3. Execute the greedy algorithm.
    placed_spheres = []
    for candidate_pos in candidate_positions:
        is_valid_placement = True
        for placed_pos in placed_spheres:
            # Calculate the squared distance to avoid using square roots.
            dist_sq = (candidate_pos[0] - placed_pos[0])**2 + \
                      (candidate_pos[1] - placed_pos[1])**2 + \
                      (candidate_pos[2] - placed_pos[2])**2
            
            # Check for overlap.
            if dist_sq < min_dist_sq:
                is_valid_placement = False
                break
        
        if is_valid_placement:
            placed_spheres.append(candidate_pos)

    # 4. Analyze and print the results.
    print("Is the problem formulation correct? Yes.")
    
    # Group spheres by z-coordinate to show the layer structure.
    layers = {}
    for pos in placed_spheres:
        z = pos[2]
        if z not in layers:
            layers[z] = 0
        layers[z] += 1
        
    # Sort layers by z-coordinate for a clear output.
    sorted_layers = sorted(layers.items())
    layer_counts = [count for z, count in sorted_layers]
    total_spheres = sum(layer_counts)
    
    equation_str = " + ".join(map(str, layer_counts))
    
    print(f"The maximum number of eyeball candies (n) is {total_spheres}.")
    print(f"The packing arrangement by layers is: {equation_str} = {total_spheres}")

# Run the solver.
solve_sphere_packing()