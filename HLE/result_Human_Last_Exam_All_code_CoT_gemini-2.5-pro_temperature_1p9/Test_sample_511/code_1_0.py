import math

def solve_packing_problem():
    """
    Searches for an optimal box for storing spherical energy balls.
    """
    
    # Initial conditions
    initial_dim = 12.0
    initial_sa = 6 * (initial_dim ** 2)
    min_balls_required = 27
    ball_radius = 2.0
    ball_diameter = 4.0

    # Best solution found so far
    min_sa_found = initial_sa
    best_dims = None

    # Search space for dimensions (L, W, H), assuming L >= W >= H
    # Steps are 0.5 cm due to manufacturing precision
    # Dimensions are unlikely to be very small or extremely large.
    h_range = [h * 0.5 for h in range(int(2 * ball_diameter), int(2 * 25))] # 4.0cm to 24.5cm
    w_range = [w * 0.5 for w in range(int(2 * ball_diameter), int(2 * 25))] # 4.0cm to 24.5cm
    l_range = [l * 0.5 for l in range(int(2 * ball_diameter), int(2 * 30))] # 4.0cm to 29.5cm
    
    for h in h_range:
        for w in w_range:
            if w < h:
                continue
            for l in l_range:
                if l < w:
                    continue
                
                # Pruning the search: if volume is too small or too large, skip
                volume = l * w * h
                min_possible_volume = min_balls_required * (4/3) * math.pi * (ball_radius**3) / 0.74 # FCC packing density
                if volume < min_possible_volume:
                    continue

                sa = 2 * (l * w + l * h + w * h)
                
                # Pruning the search: if surface area is already worse than the best, skip
                if sa >= min_sa_found:
                    continue
                
                # --- Calculate max balls for the given dimensions (l, w, h) ---
                
                # Model 1: Simple Cubic Packing
                if l < ball_diameter or w < ball_diameter or h < ball_diameter:
                    n_cubic = 0
                else:
                    nx = math.floor((l - ball_diameter) / ball_diameter) + 1
                    ny = math.floor((w - ball_diameter) / ball_diameter) + 1
                    nz = math.floor((h - ball_diameter) / ball_diameter) + 1
                    n_cubic = nx * ny * nz

                # Model 2: Staggered Layer Packing (A-B-A-B...)
                # Layer Z-coords: 2, 5, 8, ... (separation of 3cm)
                if h < ball_diameter:
                    num_z_layers = 0
                else:
                    num_z_layers = math.floor((h - ball_diameter) / (math.sqrt(8)*ball_radius/ball_diameter*2)) + 1 # sqrt(8)/2*d_ball separation
                    num_z_layers = math.floor((h - ball_diameter) / 3.0) + 1


                if num_z_layers == 0:
                    n_staggered = 0
                else:
                    num_a_layers = math.ceil(num_z_layers / 2.0)
                    num_b_layers = math.floor(num_z_layers / 2.0)
                    
                    # Balls in A layer (centers start at 2cm from edge)
                    if l < ball_diameter or w < ball_diameter:
                        balls_a = 0
                    else:
                        nx_a = math.floor((l - ball_diameter) / ball_diameter) + 1
                        ny_a = math.floor((w - ball_diameter) / ball_diameter) + 1
                        balls_a = nx_a * ny_a

                    # Balls in B layer (shifted, centers start at 4cm from edge)
                    if l < ball_diameter + 2*ball_radius or w < ball_diameter + 2*ball_radius:
                        balls_b = 0
                    else:
                        nx_b = math.floor((l - ball_radius*3) / ball_diameter) + 1 # l-6
                        ny_b = math.floor((w - ball_radius*3) / ball_diameter) + 1 # w-6
                        balls_b = nx_b * ny_b if nx_b > 0 and ny_b > 0 else 0
                    
                    n_staggered = num_a_layers * balls_a + num_b_layers * balls_b
                    
                max_balls = max(n_cubic, n_staggered)
                
                if max_balls >= min_balls_required:
                    if sa < min_sa_found:
                        min_sa_found = sa
                        # Sorting dimensions for consistent output
                        best_dims = sorted([l, w, h])

    if best_dims:
        # Check if dimensions happen to be integers
        if all(dim == int(dim) for dim in best_dims):
            dims_str = [str(int(d)) for d in best_dims]
        else:
            dims_str = [str(d) for d in best_dims]

        # Check if SA is an integer
        if min_sa_found == int(min_sa_found):
            sa_str = str(int(min_sa_found))
        else:
            sa_str = str(min_sa_found)

        print(f"{dims_str[0]}:{dims_str[1]}:{dims_str[2]}:{sa_str}")
    else:
        print(0)

# Run the solver
solve_packing_problem()