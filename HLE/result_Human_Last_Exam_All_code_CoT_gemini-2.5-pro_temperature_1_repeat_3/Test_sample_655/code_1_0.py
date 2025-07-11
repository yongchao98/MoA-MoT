import math

def solve_packing_problem():
    """
    Finds the optimal container and ball packing to maximize energy.
    This function tests three promising container candidates and uses a greedy packing
    algorithm to determine the number of balls that can fit in each.
    """

    def pack_container(container_type, dims):
        """
        Performs a greedy packing of balls into a specified container.

        Args:
            container_type (str): 'sphere', 'box', or 'cylinder'.
            dims (dict): The dimensions of the container.

        Returns:
            tuple: (number of 1-cm balls, number of 2-cm balls).
        """
        placed_balls = []
        step = 0.5

        # The algorithm iterates through two phases:
        # 1. Pack large balls (radius 2.0 cm)
        # 2. Pack small balls (radius 1.0 cm)
        # For each phase, it generates a list of all possible center points on a 0.5cm grid,
        # sorts them from the container's center outwards, and places a ball if it doesn't overlap.

        for ball_radius in [2.0, 1.0]:
            candidate_points = []
            
            # Generate a grid of valid candidate center points for the current ball size
            if container_type == 'sphere':
                r_container = dims['r']
                search_r = r_container - ball_radius
                if search_r < 0: continue
                grid_max = int(search_r / step)
                for i in range(-grid_max, grid_max + 1):
                    x = i * step
                    for j in range(-grid_max, grid_max + 1):
                        y = j * step
                        for k in range(-grid_max, grid_max + 1):
                            z = k * step
                            if x**2 + y**2 + z**2 <= search_r**2 + 1e-9:
                                candidate_points.append({'x': x, 'y': y, 'z': z})
            elif container_type == 'box':
                l, w, h = dims['l'], dims['w'], dims['h']
                max_x, max_y, max_z = l/2 - ball_radius, w/2 - ball_radius, h/2 - ball_radius
                if min(max_x, max_y, max_z) < 0: continue
                grid_x, grid_y, grid_z = int(max_x / step), int(max_y / step), int(max_z / step)
                for i in range(-grid_x, grid_x + 1):
                    x = i * step
                    for j in range(-grid_y, grid_y + 1):
                        y = j * step
                        for k in range(-grid_z, grid_z + 1):
                            z = k * step
                            candidate_points.append({'x': x, 'y': y, 'z': z})
            elif container_type == 'cylinder':
                r_container, h_container = dims['r'], dims['h']
                search_r, search_h = r_container - ball_radius, h_container/2 - ball_radius
                if search_r < 0 or search_h < 0: continue
                grid_r, grid_h = int(search_r / step), int(search_h / step)
                for i in range(-grid_r, grid_r + 1):
                    x = i * step
                    for j in range(-grid_r, grid_r + 1):
                        y = j * step
                        if x**2 + y**2 <= search_r**2 + 1e-9:
                            for k in range(-grid_h, grid_h + 1):
                                z = k * step
                                candidate_points.append({'x': x, 'y': y, 'z': z})

            # Sort points from the center outwards for a more efficient packing
            candidate_points.sort(key=lambda p: p['x']**2 + p['y']**2 + p['z']**2)

            current_phase_placed_count = 0
            for p_new in candidate_points:
                can_place = True
                min_dist_sq_check = (ball_radius * 2)**2

                # Check for overlap with existing balls
                for ball in placed_balls:
                    p_old = ball['center']
                    dist_sq = (p_new['x'] - p_old['x'])**2 + (p_new['y'] - p_old['y'])**2 + (p_new['z'] - p_old['z'])**2
                    min_dist_sq = (ball_radius + ball['radius'])**2
                    if dist_sq < min_dist_sq - 1e-9:
                        can_place = False
                        break
                
                if can_place:
                    placed_balls.append({'center': p_new, 'radius': ball_radius})

        # Count the number of balls of each type
        num_large_balls = sum(1 for b in placed_balls if b['radius'] == 2.0)
        num_small_balls = sum(1 for b in placed_balls if b['radius'] == 1.0)
        
        return num_small_balls, num_large_balls

    # --- Main Execution ---
    # Define the candidate containers to test
    candidates = [
        {'type': 'sphere', 'dims': {'r': 9.0}, 'desc': 'sphere r=9.0'},
        {'type': 'box', 'dims': {'l': 13.0, 'w': 13.0, 'h': 13.0}, 'desc': 'box 13.0x13.0x13.0'},
        {'type': 'cylinder', 'dims': {'r': 7.5, 'h': 14.5}, 'desc': 'cylinder r=7.5, h=14.5'}
    ]

    best_energy = -1
    best_config = None

    for cand in candidates:
        n1, n2 = pack_container(cand['type'], cand['dims'])
        energy = 10 * n2 + 1 * n1
        
        if energy > best_energy:
            best_energy = energy
            best_config = {
                'desc': cand['desc'],
                'n1': n1,
                'n2': n2,
                'energy': energy
            }
            
    # Print the final result
    n1 = best_config['n1']
    n2 = best_config['n2']
    energy = best_config['energy']
    desc = best_config['desc']
    
    print(f"Optimal configuration found:")
    print(f"Container: {desc}")
    print(f"Number of 1-cm balls (1 MJ): {n1}")
    print(f"Number of 2-cm balls (10 MJ): {n2}")
    print(f"Energy Calculation: 1 * {n1} + 10 * {n2} = {energy} MJ")
    
    # Final answer in the required format
    print("\nFinal Answer:")
    print(f"<<<[{desc}]{n1};{n2}>>>")

solve_packing_problem()