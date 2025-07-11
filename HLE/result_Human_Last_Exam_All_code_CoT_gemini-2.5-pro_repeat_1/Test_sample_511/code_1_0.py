import math

def solve_and_print_answer():
    """
    Analyzes packing configurations to find a box with less surface area
    than the initial 12x12x12 box that holds at least 27 balls of 4cm diameter.
    The final box dimensions must be integers.
    """
    # Initial configuration parameters
    initial_dim = 12
    ball_diameter = 4.0
    
    # Calculate initial number of balls and surface area
    n_initial = math.floor(initial_dim / ball_diameter)
    N_initial = int(n_initial ** 3)
    SA_initial = 6 * (initial_dim ** 2)

    min_sa = SA_initial
    best_config = None

    # --- Search Method 1: Simple Cubic Packing ---
    # Dimensions L, W, H must be integers. L=4*n_l, W=4*n_w, H=4*n_h
    # Search through numbers of balls (n_l, n_w, n_h) in each dimension.
    search_range = range(1, N_initial + 1)
    for n_l in search_range:
        for n_w in range(1, n_l + 1):
            # To avoid excessive computation, we can derive n_h
            if n_l * n_w == 0: continue
            n_h_min = math.ceil(N_initial / (n_l * n_w))
            for n_h in range(max(1, int(n_h_min) - 1), n_w + 1):
                if n_h == 0: continue
                num_balls = n_l * n_w * n_h
                if num_balls >= N_initial:
                    l, w, h = 4 * n_l, 4 * n_w, 4 * n_h
                    sa = 2 * (l*w + l*h + w*h)
                    if sa < min_sa:
                        min_sa = sa
                        best_config = {"l": l, "w": w, "h": h, "sa": sa}

    # --- Search Method 2: ABAB... Staggered Packing (HCP-like) ---
    # Dimensions must be integers. L=4*n_x, W=4*n_y, H=3*n_z+1
    search_range_xy = range(1, 10)
    search_range_z = range(1, 15)
    for n_x in search_range_xy:
        for n_y in range(1, n_x + 1):
            for n_z in search_range_z:
                num_balls_A = math.ceil(n_z / 2.0) * n_x * n_y
                num_balls_B = math.floor(n_z / 2.0) * (n_x - 1) * (n_y - 1) if n_x > 1 and n_y > 1 else 0
                num_balls = num_balls_A + num_balls_B
                if num_balls >= N_initial:
                    l = 4 * n_x
                    w = 4 * n_y
                    h = 3 * n_z + 1
                    sa = 2 * (l*w + l*h + w*h)
                    if sa < min_sa:
                        min_sa = sa
                        best_config = {"l": l, "w": w, "h": h, "sa": sa}

    # --- Final Output ---
    if best_config:
        # Sort dimensions for a canonical representation l>=w>=h
        dims = sorted([best_config['l'], best_config['w'], best_config['h']], reverse=True)
        l, w, h = dims[0], dims[1], dims[2]
        d = int(round(best_config['sa']))
        print(f"{l}:{w}:{h}:{d}")
    else:
        print("0")

solve_and_print_answer()