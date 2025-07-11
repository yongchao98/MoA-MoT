import math

def calculate_packing_in_box(box_dims, s, h):
    """
    Calculates the maximum number of balls that can be packed in a box
    with a given layered packing strategy (s, h).
    The function tries all 3 orientations of the box.
    """
    max_balls = 0
    # Sort dimensions to handle permutations by rotating the box
    dims = sorted(box_dims)
    from itertools import permutations
    
    # Check all 3 unique orientations (Lx, Ly, Lz), (Ly, Lz, Lx), (Lz, Lx, Ly)
    for p in set(permutations(dims)):
        L, W, H = p
        
        # Calculate how many layers fit along the height H
        # Ball diameter is 4. Center-to-wall distance is 2.
        # Available height for centers is H - 4.
        # Layer-to-layer center distance is h.
        if H - 4 < 0: continue
        num_layers = math.floor((H - 4) / h) + 1

        num_a_layers = math.ceil(num_layers / 2.0)
        num_b_layers = math.floor(num_layers / 2.0)

        # Calculate balls in A layer (the base grid)
        if L - 4 < 0 or W - 4 < 0: continue
        balls_in_a_L = math.floor((L - 4) / s) + 1
        balls_in_a_W = math.floor((W - 4) / s) + 1
        balls_in_a = balls_in_a_L * balls_in_a_W
        
        # Calculate balls in B layer (in the hollows)
        balls_in_b_L = max(0, balls_in_a_L - 1)
        balls_in_b_W = max(0, balls_in_a_W - 1)
        balls_in_b = balls_in_b_L * balls_in_b_W
        
        total_balls = num_a_layers * balls_in_a + num_b_layers * balls_in_b
        if total_balls > max_balls:
            max_balls = total_balls
            
    return max_balls

def find_best_container():
    """
    Searches for a cuboid container more efficient than the original 12x12x12 box.
    """
    original_sa = 6 * (12**2)
    min_balls = 27
    ball_diameter = 4.0

    best_solution = None
    min_sa_found = original_sa

    # Define valid staggered packing parameters (s, h)
    # s = horizontal spacing, h = vertical layer spacing
    # Constraint: (s/2)^2 + (s/2)^2 + h^2 >= D^2 => s^2/2 + h^2 >= 16
    packing_params = [
        (4.0, 3.0), # s^2/2+h^2 = 8+9=17>=16
        (5.0, 2.0), # s^2/2+h^2 = 12.5+4=16.5>=16
    ]

    # Search a range of dimensions (multiples of 0.5)
    # The range is chosen to be around the original 12x12x12 cube
    start_dim = 10.0
    end_dim = 14.0
    step = 0.5
    
    dim_range = [start_dim + i * step for i in range(int((end_dim - start_dim) / step) + 1)]

    for l in dim_range:
        for w in dim_range:
            if w < l: continue # Avoid duplicate checks (e.g., 10x11x12 vs 11x10x12)
            for h_dim in dim_range:
                if h_dim < w: continue

                sa = 2 * (l*w + l*h_dim + w*h_dim)

                if sa < min_sa_found:
                    max_balls_for_box = 0
                    for s_param, h_param in packing_params:
                        balls = calculate_packing_in_box([l, w, h_dim], s_param, h_param)
                        if balls > max_balls_for_box:
                            max_balls_for_box = balls

                    if max_balls_for_box >= min_balls:
                        min_sa_found = sa
                        best_solution = (sa, [l, w, h_dim])

    if best_solution:
        sa, dims = best_solution
        print(f"{sa:.1f}[box {dims[0]:.1f}x{dims[1]:.1f}x{dims[2]:.1f}]")
    else:
        print(0)

if __name__ == '__main__':
    find_best_container()
