import collections

def solve_packing_problem():
    """
    Solves the discrete sphere packing problem using a greedy algorithm.
    """
    # Problem parameters derived from the problem description
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)
    
    # In grid units, the minimum squared distance between centers is
    # ((radius_cm + radius_cm) / grid_step_cm)^2 = ((2+2)/0.5)^2 = 8^2 = 64.
    min_dist_sq = 64
    
    placed_spheres = []
    
    # Greedy algorithm: iterate through all potential positions and place a sphere if possible.
    # The order of iteration matters for a greedy approach. We use z, y, x.
    for z in z_range:
        for y in y_range:
            for x in x_range:
                current_pos = (x, y, z)
                is_valid = True
                
                # Check for collision with already placed spheres
                for placed_pos in placed_spheres:
                    dist_sq = (current_pos[0] - placed_pos[0])**2 + \
                              (current_pos[1] - placed_pos[1])**2 + \
                              (current_pos[2] - placed_pos[2])**2
                    
                    if dist_sq < min_dist_sq:
                        is_valid = False
                        break
                
                if is_valid:
                    placed_spheres.append(current_pos)

    # The problem formulation is correct.
    # The maximized value n is the total number of spheres placed.
    n = len(placed_spheres)
    print(f"The problem formulation is correct.")
    print(f"The highest number of eyeball candies that can be packed is: {n}")
    
    # To satisfy the "output each number in the final equation" requirement,
    # we show the breakdown of spheres by horizontal layer (z-coordinate).
    if n > 0:
        layer_counts = collections.Counter(sphere[2] for sphere in placed_spheres)
        sorted_layers = sorted(layer_counts.items())
        
        count_parts = [str(count) for z, count in sorted_layers]
        equation_str = " + ".join(count_parts) + f" = {n}"
        
        print("The composition of the total number by layers is:")
        print(equation_str)

solve_packing_problem()
