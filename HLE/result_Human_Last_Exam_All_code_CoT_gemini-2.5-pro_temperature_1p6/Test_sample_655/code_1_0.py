import math

def solve_packing_problem():
    """
    Main function to solve the packing problem.
    It iterates through a list of candidate containers and uses a greedy
    packing simulation to find the one that yields the maximum energy.
    """
    best_config = {'energy': -1}

    def update_best(container_type, dims_str, n1, n2):
        """Updates the best configuration found so far."""
        nonlocal best_config
        energy = 1 * n1 + 10 * n2
        if energy > best_config['energy']:
            best_config = {
                'energy': energy,
                'type': container_type,
                'dims_str': dims_str,
                'n1': n1,
                'n2': n2,
            }

    def do_pack(is_inside_func, bounds_u):
        """
        Performs a greedy packing simulation for a given container shape.
        All coordinates are in integer units of 0.5 cm.
        The function first packs larger, more valuable balls, then fills gaps with smaller balls.
        """
        placed_balls = []  # List of ((x, y, z), radius_u)

        # --- Pass 1: Pack 2-cm balls (radius = 4 units) ---
        r_u_large = 4
        min_x, max_x, min_y, max_y, min_z, max_z = bounds_u
        
        # A simple raster scan iteration is used for placing balls.
        for z in range(max_z - r_u_large, min_z + r_u_large - 1, -2):
            for y in range(max_y - r_u_large, min_y + r_u_large - 1, -2):
                for x in range(max_x - r_u_large, min_x + r_u_large - 1, -2):
                    center_u = (x, y, z)

                    if is_inside_func(center_u, r_u_large):
                        can_place = True
                        for p_center_u, p_r_u in placed_balls:
                            dist_sq = sum((c1 - c2) ** 2 for c1, c2 in zip(center_u, p_center_u))
                            if dist_sq < (r_u_large + p_r_u) ** 2:
                                can_place = False
                                break
                        if can_place:
                            placed_balls.append((center_u, r_u_large))
        n2 = len(placed_balls)

        # --- Pass 2: Pack 1-cm balls (radius = 2 units) ---
        r_u_small = 2
        for z in range(max_z - r_u_small, min_z + r_u_small - 1, -1):
            for y in range(max_y - r_u_small, min_y + r_u_small - 1, -1):
                for x in range(max_x - r_u_small, min_x + r_u_small - 1, -1):
                    center_u = (x, y, z)
                    if is_inside_func(center_u, r_u_small):
                        can_place = True
                        for p_center_u, p_r_u in placed_balls:
                            dist_sq = sum((c1 - c2) ** 2 for c1, c2 in zip(center_u, p_center_u))
                            if dist_sq < (r_u_small + p_r_u) ** 2:
                                can_place = False
                                break
                        if can_place:
                            placed_balls.append((center_u, r_u_small))
        n1 = len(placed_balls) - n2
        return n1, n2

    # A curated list of promising candidates to test.
    candidates = [
        # Spheres (theoretically most volume-efficient for a given surface area)
        {'type': 'sphere', 'R': 9.0},
        {'type': 'sphere', 'R': 8.5},
        # Cylinders (near-optimal where height is close to diameter)
        {'type': 'cylinder', 'r': 8.0, 'h': 12.0},  # A=1005 cm^2
        {'type': 'cylinder', 'r': 7.5, 'h': 13.0},  # A=964 cm^2
        # Boxes (cubes are the most volume-efficient boxes)
        {'type': 'box', 'L': 13.0, 'W': 13.0, 'H': 13.0},  # A=1014 cm^2
        {'type': 'box', 'L': 13.5, 'W': 13.0, 'H': 12.5},  # A=1013.5 cm^2
    ]

    for cand in candidates:
        if cand['type'] == 'sphere':
            R = cand['R']
            if 4 * math.pi * R**2 > 1050: continue
            R_u = int(R / 0.5)
            def is_inside_sphere(center_u, r_u):
                dist_sq_from_origin = center_u[0]**2 + center_u[1]**2 + center_u[2]**2
                return dist_sq_from_origin <= (R_u - r_u)**2
            bounds_u = (-R_u, R_u, -R_u, R_u, -R_u, R_u)
            n1, n2 = do_pack(is_inside_sphere, bounds_u)
            dims_str = f"sphere r={R:.1f}"
            update_best('sphere', dims_str, n1, n2)

        elif cand['type'] == 'cylinder':
            r, h = cand['r'], cand['h']
            if 2 * math.pi * r * (r + h) > 1050: continue
            r_u, h_u = int(r / 0.5), int(h / 0.5)
            def is_inside_cylinder(center_u, r_u):
                cx, cy, cz = center_u
                if not (r_u <= cz <= h_u - r_u): return False
                return cx**2 + cy**2 <= (r_u - r_u)**2
            bounds_u = (-r_u, r_u, -r_u, r_u, 0, h_u)
            n1, n2 = do_pack(is_inside_cylinder, bounds_u)
            dims_str = f"cylinder r={r:.1f}, h={h:.1f}"
            update_best('cylinder', dims_str, n1, n2)

        elif cand['type'] == 'box':
            L, W, H = cand['L'], cand['W'], cand['H']
            if 2 * (L*W + L*H + W*H) > 1050: continue
            L_u, W_u, H_u = int(L/0.5), int(W/0.5), int(H/0.5)
            # The check is handled by the iteration bounds
            def is_inside_box(center_u, r_u): return True
            bounds_u = (0, L_u, 0, W_u, 0, H_u)
            n1, n2 = do_pack(is_inside_box, bounds_u)
            dims_str = f"box {L:.1f}x{W:.1f}x{H:.1f}"
            update_best('box', dims_str, n1, n2)

    # Print the final result
    b = best_config
    final_answer_string = f"{b['dims_str']};{b['n1']};{b['n2']}"
    energy_calculation_string = f"Total energy: 1 * {b['n1']} + 10 * {b['n2']} = {b['energy']}"
    
    print(energy_calculation_string)
    print(f"<<<{final_answer_string}>>>")

solve_packing_problem()