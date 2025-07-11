import math

def solve_packing_problem():
    """
    Main function to find the optimal container and ball packing configuration.
    It iterates through shapes, dimensions, and uses a greedy packing algorithm.
    """
    
    # --- Constants ---
    MAX_SA = 1050.0
    STEP = 0.5
    R_LARGE, E_LARGE = 2.0, 10
    R_SMALL, E_SMALL = 1.0, 1
    MIN_DIM_FOR_LARGE_BALL = 4.0

    # --- State for tracking the best result ---
    best_result = {
        "energy": -1,
        "config_str": "",
        "n_small": 0,
        "n_large": 0
    }

    def distance_sq(p1, p2):
        """Calculates the squared distance between two 3D points."""
        return (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2

    def get_valid_centers(dims, shape, radius):
        """Generates a list of valid center coordinates within a container for a given ball radius."""
        centers = []
        if shape == 'box':
            L, W, H = dims
            # Valid centers must be `radius` away from each wall
            x_min_i, x_max_i = math.ceil(radius / STEP), math.floor((L - radius) / STEP)
            y_min_i, y_max_i = math.ceil(radius / STEP), math.floor((W - radius) / STEP)
            z_min_i, z_max_i = math.ceil(radius / STEP), math.floor((H - radius) / STEP)
            
            for i in range(x_min_i, x_max_i + 1):
                for j in range(y_min_i, y_max_i + 1):
                    for k in range(z_min_i, z_max_i + 1):
                        centers.append((i * STEP, j * STEP, k * STEP))
                        
        elif shape == 'cylinder':
            R, H = dims
            r_boundary = R - radius
            z_min_i, z_max_i = math.ceil(radius / STEP), math.floor((H - radius) / STEP)
            xy_max_i = math.floor(r_boundary / STEP)
            
            # To make it symmetrical and centered for easier calculations
            # We iterate in a bounding box and check the circular constraint
            for k in range(z_min_i, z_max_i + 1):
                z = k * STEP
                for i in range(-xy_max_i, xy_max_i + 1):
                    x = i * STEP
                    y_max_coord = math.sqrt(max(0, r_boundary**2 - x**2))
                    y_max_i = math.floor(y_max_coord / STEP)
                    for j in range(-y_max_i, y_max_i + 1):
                        y = j * STEP
                        # The box is defined from z=0 to H, so we shift center z
                        centers.append((x, y, z))
                        
        else: # sphere
            R = dims[0]
            r_boundary = R - radius
            r_boundary_sq = r_boundary**2
            max_i = math.floor(r_boundary / STEP)
            
            for k in range(-max_i, max_i + 1):
                z = k * STEP
                for i in range(-max_i, max_i + 1):
                    x = i * STEP
                    if x**2 + z**2 <= r_boundary_sq:
                        y_max_coord = math.sqrt(max(0, r_boundary_sq - x**2 - z**2))
                        y_max_i = math.floor(y_max_coord / STEP)
                        for j in range(-y_max_i, y_max_i + 1):
                            y = j * STEP
                            centers.append((x, y, z))
        return centers

    def pack_container(dims, shape):
        """Performs the greedy packing and updates the best result if a better one is found."""
        
        placed_balls = []
        
        # --- Part 1: Place Large Balls (r=2) ---
        centers_large = get_valid_centers(dims, shape, R_LARGE)
        centers_large.sort(key=lambda p: (p[2], p[1], p[0])) # Deterministic packing order

        for center in centers_large:
            can_place = True
            min_dist_sq_ll = (R_LARGE + R_LARGE)**2
            for p_center, _ in placed_balls:
                if distance_sq(center, p_center) < min_dist_sq_ll:
                    can_place = False
                    break
            if can_place:
                placed_balls.append((center, R_LARGE))

        num_large = len(placed_balls)
        
        # --- Part 2: Place Small Balls (r=1) ---
        centers_small = get_valid_centers(dims, shape, R_SMALL)
        centers_small.sort(key=lambda p: (p[2], p[1], p[0]))

        for center in centers_small:
            can_place = True
            for p_center, p_radius in placed_balls:
                min_dist_sq = (R_SMALL + p_radius)**2
                if distance_sq(center, p_center) < min_dist_sq:
                    can_place = False
                    break
            if can_place:
                placed_balls.append((center, R_SMALL))
                
        num_small = len(placed_balls) - num_large
        energy = num_large * E_LARGE + num_small * E_SMALL
        
        if energy > best_result["energy"]:
            best_result["energy"] = energy
            best_result["n_small"] = num_small
            best_result["n_large"] = num_large
            
            # Format dimension string, removing .0 for integers
            if shape == 'box':
                L, W, H = (int(d) if d == int(d) else d for d in dims)
                best_result["config_str"] = f"box {L}x{W}x{H}"
            elif shape == 'cylinder':
                R, H = (int(d) if d == int(d) else d for d in dims)
                best_result["config_str"] = f"cylinder r={R}, h={H}"
            else: # sphere
                R = dims[0]
                R = int(R) if R == int(R) else R
                best_result["config_str"] = f"sphere r={R}"

    # --- Iterate through Box Dimensions ---
    max_dim_cap = 20.0
    h_i_max = math.floor(min(13.5, max_dim_cap) / STEP)
    for h_i in range(math.ceil(MIN_DIM_FOR_LARGE_BALL / STEP), h_i_max + 1):
        H = h_i * STEP
        w_i_max = math.floor(min(-H + math.sqrt(H**2 + 2 * MAX_SA), max_dim_cap) / STEP)
        for w_i in range(h_i, w_i_max + 1):
            W = w_i * STEP
            if (2*W + 2*H) == 0: continue
            max_l = (MAX_SA - 2*W*H) / (2*W + 2*H)
            l_i_max = math.floor(min(max_l, max_dim_cap) / STEP)
            for l_i in range(w_i, l_i_max + 1):
                L = l_i * STEP
                if 2 * (L*W + L*H + W*H) > MAX_SA: continue
                pack_container((L, W, H), 'box')

    # --- Iterate through Cylinder Dimensions ---
    r_i_max = math.floor(math.sqrt(MAX_SA / (2 * math.pi)) / STEP)
    for r_i in range(math.ceil(R_LARGE / STEP), r_i_max + 1):
        R = r_i * STEP
        if (2 * math.pi * R) == 0: continue
        max_h = (MAX_SA / (2 * math.pi * R)) - R
        if max_h < MIN_DIM_FOR_LARGE_BALL: continue
        for h_i in range(math.ceil(MIN_DIM_FOR_LARGE_BALL / STEP), math.floor(max_h / STEP) + 1):
            H = h_i * STEP
            pack_container((R, H), 'cylinder')
            
    # --- Iterate through Sphere Dimensions ---
    r_i_max = math.floor(math.sqrt(MAX_SA / (4 * math.pi)) / STEP)
    for r_i in range(math.ceil(R_LARGE / STEP), r_i_max + 1):
        R = r_i * STEP
        pack_container((R,), 'sphere')
    
    # --- Print Final Answer ---
    C = best_result["config_str"]
    a = best_result["n_small"]
    b = best_result["n_large"]
    
    total_energy = best_result["energy"]
    
    print(f"[{C}]{a};{b}")
    print(f"Final Equation: {total_energy} = 1 * {a} + 10 * {b}")


solve_packing_problem()
<<<[sphere r=9]16;32>>>