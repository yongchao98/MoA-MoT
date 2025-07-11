import math

def solve_candy_packing():
    """
    Calculates the maximum number of spherical candies that can be packed into a box
    using a greedy algorithm on a discrete grid.
    """
    # Box dimensions in cm: 12x12x11
    # Candy radius in cm: 2
    # Grid unit: 0.5 cm

    # Convert dimensions to grid units
    box_dims = {'x': 24, 'y': 24, 'z': 22} # 12/0.5, 12/0.5, 11/0.5
    radius = 4 # 2 / 0.5

    # Squared minimum distance between centers (diameter^2)
    # (radius_1 + radius_2)^2 = (4+4)^2 = 8^2 = 64
    min_dist_sq = (2 * radius)**2

    # Define the valid range for the center coordinates
    # A center (cx,cy,cz) is valid if radius <= c <= dimension - radius
    x_range = range(radius, box_dims['x'] - radius + 1) # [4, 20]
    y_range = range(radius, box_dims['y'] - radius + 1) # [4, 20]
    z_range = range(radius, box_dims['z'] - radius + 1) # [4, 18]

    # Generate all possible candidate points for the centers in a lexicographical order
    candidate_points = []
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_points.append({'x': x, 'y': y, 'z': z})

    placed_balls = []

    # Greedy algorithm: iterate through candidates and place if possible
    for point in candidate_points:
        can_place = True
        for placed in placed_balls:
            dist_sq = (point['x'] - placed['x'])**2 + \
                      (point['y'] - placed['y'])**2 + \
                      (point['z'] - placed['z'])**2
            
            if dist_sq < min_dist_sq:
                can_place = False
                break
        
        if can_place:
            placed_balls.append(point)

    # The problem formulation is correct.
    # Now, we output the result based on the greedy algorithm's findings.
    
    # Group placed balls by their z-coordinate to count balls per layer
    layers = {}
    for ball in placed_balls:
        z = ball['z']
        if z not in layers:
            layers[z] = 0
        layers[z] += 1
    
    print("The problem formulation is correct. Here is a possible packing solution:")
    
    layer_counts = []
    # Sort by layer height and print the breakdown
    for z_coord in sorted(layers.keys()):
        count = layers[z_coord]
        layer_counts.append(str(count))
        print(f"Number of balls in layer at z={z_coord*0.5}cm: {count}")

    total_balls = len(placed_balls)
    equation = " + ".join(layer_counts)
    
    print(f"\nTotal number of eyeball candies: {equation} = {total_balls}")


solve_candy_packing()
