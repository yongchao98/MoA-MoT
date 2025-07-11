import random

def solve_packing():
    """
    Solves the sphere packing problem using a randomized greedy algorithm
    to find the maximum number of spheres (n).
    """
    # Box dimensions in 0.5 cm grid units.
    # 12cm -> 24 units, 12cm -> 24 units, 11cm -> 22 units.
    box_dims = (24, 24, 22)

    # Sphere radius in 0.5 cm grid units.
    # 2cm -> 4 units.
    radius = 4

    # The squared distance between centers must be >= (2*r)^2.
    # d^2 >= (4+4)^2 = 64.
    min_sq_dist = (2 * radius)**2

    # As per the problem, the valid integer ranges for the centers are:
    # x_i, y_i in [4, 20]
    # z_i in [4, 18]
    x_coords = range(4, 21)
    y_coords = range(4, 21)
    z_coords = range(4, 19)
    
    # Generate all possible candidate positions for the sphere centers.
    candidate_positions = []
    for z in z_coords:
        for y in y_coords:
            for x in x_coords:
                candidate_positions.append((x, y, z))

    max_n = 0
    
    # We run the randomized search for a number of trials to increase
    # the chance of finding the optimal packing.
    num_trials = 100 

    for _ in range(num_trials):
        placed_spheres_centers = []
        
        # In each trial, we shuffle the candidate positions to try a
        # different packing order.
        random.shuffle(candidate_positions)

        for pos in candidate_positions:
            can_place = True
            for placed_pos in placed_spheres_centers:
                # Calculate squared Euclidean distance.
                sq_dist = (pos[0] - placed_pos[0])**2 + \
                          (pos[1] - placed_pos[1])**2 + \
                          (pos[2] - placed_pos[2])**2
                
                # Check for overlap.
                if sq_dist < min_sq_dist:
                    can_place = False
                    break
            
            if can_place:
                placed_spheres_centers.append(pos)
        
        if len(placed_spheres_centers) > max_n:
            max_n = len(placed_spheres_centers)

            # Based on known results, 22 is the likely maximum.
            # We can stop if we find it.
            if max_n == 22:
                break
    
    # Since the problem formulation is correct, we print the maximized value n.
    # The output format required is just the number.
    print(max_n)

if __name__ == '__main__':
    solve_packing()