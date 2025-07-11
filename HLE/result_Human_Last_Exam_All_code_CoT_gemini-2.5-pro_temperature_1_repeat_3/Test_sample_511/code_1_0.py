import math

def calculate_packing_capacity(length, width, height):
    """
    Calculates the number of spheres (radius=2) that can be packed into a box.
    This function implements a dense packing model (a variant of HCP/FCC)
    with centers on a 0.5cm grid.
    It assumes the layers are stacked along the 'height' dimension.
    
    - Layer A consists of a square grid of balls.
    - Layer B is another square grid, shifted to fit in the hollows of layer A.
    - The layers are stacked in an alternating pattern (A, B, A, B, ...).
    """
    if any(d < 4 for d in [length, width, height]):
        return 0

    # Layer A centers: (2+4i, 2+4j, z_k)
    # Number of balls in a full "A" layer
    # The first center coordinate is 2cm from the edge.
    # The space available for centers is (dimension - 2*radius).
    # The distance between centers is 4cm (diameter).
    nx_A = math.floor((length - 4) / 4) + 1
    ny_A = math.floor((width - 4) / 4) + 1
    balls_in_layer_A = nx_A * ny_A
    
    # Layer B centers: (4+4i, 4+4j, z_k + 3)
    # This offset is chosen to maximize density while keeping centers on the 0.5cm grid
    # and ensuring distance between balls in adjacent layers is >= 4cm.
    # sqrt((4-2)^2 + (4-2)^2 + 3^2) = sqrt(4+4+9) = sqrt(17) >= 4
    if length < 6 or width < 6:
        balls_in_layer_B = 0
    else:
        nx_B = math.floor((length - 6) / 4) + 1
        ny_B = math.floor((width - 6) / 4) + 1
        balls_in_layer_B = nx_B * ny_B
        
    # Number of layers along the height dimension. Layer separation is 3cm.
    num_layers = math.floor((height - 4) / 3) + 1
    if num_layers <= 0:
        return 0
        
    # The number of A and B layers depends on the total number of layers.
    num_A_layers = math.ceil(num_layers / 2)
    num_B_layers = math.floor(num_layers / 2)
    
    total_balls = num_A_layers * balls_in_layer_A + num_B_layers * balls_in_layer_B
    return total_balls

def find_optimal_box():
    """
    Searches for a box with integer dimensions that has a smaller surface area
    than the initial 12x12x12 box but holds at least the same number of balls.
    """
    initial_balls = 27
    initial_surface_area = 864

    best_dims = None
    min_surface_area = initial_surface_area

    # Search range for integer dimensions (in cm).
    # To have S < 864, at least one dimension must be < 12.
    # We loop L >= W >= H to avoid redundant permutations.
    # From 6H^2 <= S < 864, we get H < 12. So H is in [4, 11].
    for h in range(4, 12):
        # From 2(W*H + W*H + H*H) <= 2(L*W + L*H + W*H) < 864, we get 2(2W*H + H^2) < 864
        # W < (432 - H^2)/(2H)
        w_limit = math.ceil((432 - h**2) / (2 * h))
        for w in range(h, w_limit):
            # Pruning the search: if the area of the two largest faces already exceeds the minimum, stop.
            if 2 * w * h >= min_surface_area:
                break
            # From 2(L*W + L*H + W*H) < 864, we get L < (432 - W*H)/(W+H)
            l_limit = math.ceil((432 - w * h) / (w + h))
            for l in range(w, l_limit):
                current_surface_area = 2 * (l * w + l * h + w * h)

                if current_surface_area >= min_surface_area:
                    continue

                # We must check all 3 possible orientations for packing, as the
                # packing function is dependent on which dimension is 'height'.
                # The dimensions are (l, w, h).
                balls_p1 = calculate_packing_capacity(l, w, h) # stack along h
                balls_p2 = calculate_packing_capacity(l, h, w) # stack along w
                balls_p3 = calculate_packing_capacity(w, h, l) # stack along l
                
                max_balls = max(balls_p1, balls_p2, balls_p3)

                if max_balls >= initial_balls:
                    # Found a better solution
                    min_surface_area = current_surface_area
                    # Store dimensions sorted from largest to smallest
                    best_dims = tuple(sorted((l, w, h), reverse=True))

    if best_dims:
        l, w, h = best_dims
        s = int(min_surface_area)
        print(f"{l}:{w}:{h}:{s}")
    else:
        print(0)

find_optimal_box()