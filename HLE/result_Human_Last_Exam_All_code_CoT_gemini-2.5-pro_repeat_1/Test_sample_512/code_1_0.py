import math

def calculate_packing(L, W, H):
    """
    Calculates the maximum number of balls that can fit in a L x W x H box
    using a staggered layer packing (A-B layers).
    The packing is tested along all three axes.
    """
    
    # Radius of ball is 2, so center must be >= 2 from any wall.
    # Min distance between centers is 4.
    # Min vertical distance between A and B layer centers is 3.

    max_balls = 0
    
    # Try layering along each axis (H, W, L)
    for dim_base1, dim_base2, dim_height in [(L, W, H), (L, H, W), (W, H, L)]:
        
        if dim_height < 4:
            continue

        balls_A = math.floor(dim_base1 / 4) * math.floor(dim_base2 / 4)
        balls_B = math.floor((dim_base1 - 2) / 4) * math.floor((dim_base2 - 2) / 4)
        
        if balls_A == 0:
            continue

        z_center_max = dim_height - 2.0
        
        current_z = 2.0
        total_balls = 0
        is_layer_A = True
        
        while current_z <= z_center_max:
            if is_layer_A:
                total_balls += balls_A
            else:
                total_balls += balls_B
            
            is_layer_A = not is_layer_A
            current_z += 3.0 # Staggered layers can be 3cm apart

        if total_balls > max_balls:
            max_balls = total_balls
            
    return max_balls

def find_best_container():
    """
    Searches for the container with the minimum surface area that can hold at least 27 balls.
    """
    initial_area = 864
    min_balls = 27
    
    best_solution = None
    min_area = initial_area
    
    # Define a reasonable search range for dimensions (in multiples of 0.5 cm)
    # A side length > 12 would likely increase surface area unless other sides are very small.
    l_range = [i * 0.5 for i in range(8, 28)] # L from 4.0 to 14.0
    w_range = [i * 0.5 for i in range(8, 28)] # W from 4.0 to 14.0
    h_range = [i * 0.5 for i in range(8, 28)] # H from 4.0 to 14.0

    for l in l_range:
        for w in w_range:
            # To avoid duplicate checks like 10x11x12 and 11x10x12, ensure w>=l
            if w < l:
                continue
            for h in h_range:
                # Ensure h>=w
                if h < w:
                    continue

                area = 2 * (l*w + l*h + w*h)
                
                if area >= min_area:
                    continue
                
                num_balls = calculate_packing(l, w, h)
                
                if num_balls >= min_balls:
                    min_area = area
                    best_solution = {
                        "area": area,
                        "dims": (l, w, h),
                        "balls": num_balls,
                        "type": "box"
                    }

    if best_solution:
        dims = best_solution["dims"]
        # Format dimensions to remove trailing .0
        l_str = f"{dims[0]:.1f}".replace('.0', '')
        w_str = f"{dims[1]:.1f}".replace('.0', '')
        h_str = f"{dims[2]:.1f}".replace('.0', '')
        
        area_str = f"{best_solution['area']:.1f}".replace('.0', '')

        print(f"{area_str}[box {l_str}x{w_str}x{h_str}]")
    else:
        print("0")

find_best_container()