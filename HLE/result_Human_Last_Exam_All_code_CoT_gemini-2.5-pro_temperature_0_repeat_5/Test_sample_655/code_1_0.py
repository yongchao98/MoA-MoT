import math

def solve():
    """
    Finds the optimal container and ball packing to maximize energy.
    """
    best_result = {"energy": -1, "n1": 0, "n2": 0, "desc": ""}

    def pack_container(shape, dims):
        # This function simulates the packing process for a given container.
        # It uses integer coordinates where 1 unit = 0.5 cm.
        # Radii in integer units: r_1cm -> 2, r_2cm -> 4.
        # Non-overlap conditions are checked using squared distances to avoid sqrt.
        # dist_sq_int >= 4 * (r1+r2)^2
        # r1+r1: 4*(1+1)^2 = 16
        # r1+r2: 4*(1+2)^2 = 36
        # r2+r2: 4*(2+2)^2 = 64
        
        memo_inside = {}
        dims_int = {k: int(v * 2) for k, v in dims.items()}

        def is_inside(p_int, r_ball_int, shape, dims_int):
            key = (p_int, r_ball_int, shape, tuple(sorted(dims_int.items())))
            if key in memo_inside:
                return memo_inside[key]

            ix, iy, iz = p_int
            res = False
            if shape == 'sphere':
                # Ball is inside if dist from origin + ball_radius <= container_radius
                res = ix**2 + iy**2 + iz**2 <= (dims_int['r'] - r_ball_int)**2
            elif shape == 'box':
                # Using a centered coordinate system [-L/2, L/2]
                res = (abs(ix) + r_ball_int <= dims_int['L'] / 2 and
                       abs(iy) + r_ball_int <= dims_int['W'] / 2 and
                       abs(iz) + r_ball_int <= dims_int['H'] / 2)
            elif shape == 'cylinder':
                # Using a centered coordinate system, axis on z
                res = (ix**2 + iy**2 <= (dims_int['r'] - r_ball_int)**2 and
                       abs(iz) + r_ball_int <= dims_int['h'] / 2)
            memo_inside[key] = res
            return res

        # Define bounding box for grid search based on centered coordinates
        if shape == 'sphere':
            max_dim = dims_int['r']
            min_x, max_x = -max_dim, max_dim
            min_y, max_y = -max_dim, max_dim
            min_z, max_z = -max_dim, max_dim
        elif shape == 'box':
            min_x, max_x = -dims_int['L'] // 2, dims_int['L'] // 2
            min_y, max_y = -dims_int['W'] // 2, dims_int['W'] // 2
            min_z, max_z = -dims_int['H'] // 2, dims_int['H'] // 2
        elif shape == 'cylinder':
            min_x, max_x = -dims_int['r'], dims_int['r']
            min_y, max_y = -dims_int['r'], dims_int['r']
            min_z, max_z = -dims_int['h'] // 2, dims_int['h'] // 2

        all_points = []
        for ix in range(min_x, max_x + 1):
            for iy in range(min_y, max_y + 1):
                for iz in range(min_z, max_z + 1):
                    all_points.append((ix, iy, iz))

        # --- Pass 1: Place 2-cm radius balls (r_int = 4) ---
        r2_int = 4
        valid_centers_r2 = [p for p in all_points if is_inside(p, r2_int, shape, dims_int)]
        valid_centers_r2.sort(key=lambda p: (p[0]**2 + p[1]**2 + p[2]**2, p[0], p[1], p[2]))

        placed_balls = []
        for p_new in valid_centers_r2:
            can_place = True
            for p_placed, r_placed in placed_balls:
                dist_sq = (p_new[0] - p_placed[0])**2 + (p_new[1] - p_placed[1])**2 + (p_new[2] - p_placed[2])**2
                if dist_sq < 64:
                    can_place = False
                    break
            if can_place:
                placed_balls.append((p_new, r2_int))
        
        n2 = len(placed_balls)

        # --- Pass 2: Place 1-cm radius balls (r_int = 2) ---
        r1_int = 2
        valid_centers_r1 = [p for p in all_points if is_inside(p, r1_int, shape, dims_int)]
        valid_centers_r1.sort(key=lambda p: (p[0]**2 + p[1]**2 + p[2]**2, p[0], p[1], p[2]))
        
        n1_added = 0
        for p_new in valid_centers_r1:
            can_place = True
            for p_placed, r_placed_int in placed_balls:
                dist_sq = (p_new[0] - p_placed[0])**2 + (p_new[1] - p_placed[1])**2 + (p_new[2] - p_placed[2])**2
                req_dist_sq = 16 if r_placed_int == r1_int else 36
                if dist_sq < req_dist_sq:
                    can_place = False
                    break
            if can_place:
                placed_balls.append((p_new, r1_int))
                n1_added += 1
        
        n1 = n1_added
        return 10 * n2 + 1 * n1, n1, n2

    # Define promising container configurations to test
    configs = []
    max_sa = 1050.0

    # Spheres (r from 4.0 to 9.0)
    for r_int in range(8, 19):
        r = r_int / 2.0
        if 4 * math.pi * r**2 <= max_sa:
            configs.append({'shape': 'sphere', 'dims': {'r': r}, 'desc': f"sphere r={r}"})

    # Cubes (L from 4.0 to 13.0)
    for L_int in range(8, 27):
        L = L_int / 2.0
        if 6 * L**2 <= max_sa:
            configs.append({'shape': 'box', 'dims': {'L': L, 'W': L, 'H': L}, 'desc': f"box {L}x{L}x{L}"})

    # Cylinders with h=2r (r from 4.0 to 7.0)
    for r_int in range(8, 15):
        r = r_int / 2.0
        h = 2 * r
        if (2 * math.pi * r**2 + 2 * math.pi * r * h) <= max_sa:
            configs.append({'shape': 'cylinder', 'dims': {'r': r, 'h': h}, 'desc': f"cylinder r={r}, h={h}"})
            
    # A few other promising cylinders (pancakes and rods)
    other_cylinders = [{'r': 9.0, 'h': 9.5}, {'r': 8.5, 'h': 11.0}, {'r': 4.0, 'h': 38.0}]
    for d in other_cylinders:
        r, h = d['r'], d['h']
        if (2 * math.pi * r**2 + 2 * math.pi * r * h) <= max_sa:
             configs.append({'shape': 'cylinder', 'dims': d, 'desc': f"cylinder r={r}, h={h}"})

    # Run the packing simulation for all configurations
    for config in configs:
        energy, n1, n2 = pack_container(config['shape'], config['dims'])
        if energy > best_result['energy']:
            best_result['energy'] = energy
            best_result['n1'] = n1
            best_result['n2'] = n2
            best_result['desc'] = config['desc']

    final_desc = best_result['desc']
    final_n1 = best_result['n1']
    final_n2 = best_result['n2']
    
    print(f"[{final_desc}]{final_n1};{final_n2}")

solve()
<<<[box 13.0x13.0x13.0]8;27>>>