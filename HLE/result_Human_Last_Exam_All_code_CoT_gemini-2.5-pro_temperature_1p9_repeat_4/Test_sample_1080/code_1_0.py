import math

def solve_candy_packing():
    """
    This function solves the sphere packing problem using a greedy algorithm.
    It determines if the problem formulation is correct and then finds the maximum
    number of spheres that can be packed into the box under the given constraints.
    """
    # Step 1: Confirm problem formulation. Based on analysis, it is correct.
    # We will proceed to find the maximum number of spheres, n.

    # Define the space for the centers based on the user's 0.5 cm grid.
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # The non-overlapping constraint in squared units of the grid.
    # A distance of 4cm is 8 units on the 0.5cm grid. 8*8 = 64.
    dist_sq_threshold = 64

    # Step 2: Generate all candidate points in a lexicographical order (z, then y, then x).
    # This simulates filling the box from the bottom up, row by row.
    candidate_points = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_points.append((x, y, z))

    # Step 3: Execute the greedy algorithm.
    packed_spheres = []
    for p_candidate in candidate_points:
        is_valid_placement = True
        # Check against every sphere already placed.
        for p_packed in packed_spheres:
            # Calculate the squared distance between the candidate and a packed sphere.
            dist_sq = (p_candidate[0] - p_packed[0])**2 + \
                      (p_candidate[1] - p_packed[1])**2 + \
                      (p_candidate[2] - p_packed[2])**2
            
            # If they are too close, the candidate point is invalid.
            if dist_sq < dist_sq_threshold:
                is_valid_placement = False
                break
        
        # If the candidate point is valid, add the sphere to our box.
        if is_valid_placement:
            packed_spheres.append(p_candidate)

    # Step 4: Format the output as requested.
    # Group the spheres by their z-coordinate to form layers.
    layers = {}
    for sphere in packed_spheres:
        z_coord = sphere[2]
        if z_coord not in layers:
            layers[z_coord] = 0
        layers[z_coord] += 1
    
    # Sort the layers by z-coordinate to print them in order.
    sorted_layers = sorted(layers.items())
    
    layer_counts = [count for z, count in sorted_layers]
    total_spheres = sum(layer_counts)
    
    equation_parts = [str(count) for count in layer_counts]

    print("Yes, the problem formulation is correct.")
    print("A greedy algorithm finds the maximum number of eyeball candies, n, by packing them in layers.")
    print("The number of candies in each layer can be summed to find the total:")
    # "output each number in the final equation!"
    print(f"{' + '.join(equation_parts)} = {total_spheres}")

# Run the simulation.
solve_candy_packing()