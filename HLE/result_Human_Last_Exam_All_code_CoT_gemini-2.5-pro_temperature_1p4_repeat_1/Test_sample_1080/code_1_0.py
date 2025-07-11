import math

def solve_packing():
    """
    Solves the sphere packing problem for the given box and sphere dimensions
    using a greedy algorithm.
    """
    # Problem parameters in integer grid units (0.5 cm)
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)
    # The radius is 2cm, which is 4 units on the 0.5cm grid.
    # The squared distance between centers must be >= (2 * radius)^2 = (2 * 4)^2 = 64.
    min_sq_dist = 64

    placed_spheres = []
    spheres_per_layer = {}

    # Iterate through each possible position, layering along the shortest (z) axis
    for z in z_range:
        count_in_layer = 0
        for y in y_range:
            for x in x_range:
                candidate_center = (x, y, z)
                can_place = True
                
                # Check for overlap with already placed spheres
                for placed_center in placed_spheres:
                    sq_dist = ( (candidate_center[0] - placed_center[0]) ** 2 +
                                (candidate_center[1] - placed_center[1]) ** 2 +
                                (candidate_center[2] - placed_center[2]) ** 2 )
                    
                    if sq_dist < min_sq_dist:
                        can_place = False
                        break
                
                if can_place:
                    placed_spheres.append(candidate_center)
                    count_in_layer += 1
        
        if count_in_layer > 0:
            spheres_per_layer[z] = count_in_layer
            
    # Prepare the final output equation
    layer_counts = list(spheres_per_layer.values())
    total_spheres = sum(layer_counts)
    
    equation = " + ".join(map(str, layer_counts))
    
    # Print each number in the final equation
    print(f"{equation} = {total_spheres}")

# Run the solver
solve_packing()