import math

def calculate_packed_balls(l, w, h):
    """
    Calculates the number of balls that can be packed in a box of size L x W x H
    using a staggered A-B-A-B packing.
    """
    if l < 4 or w < 4 or h < 4:
        return 0

    # Calculate how many layers fit in the height H
    # Layer centers are at z = 2, 5, 8, ...
    # Max center z is limited by h: z_center + radius <= h => z_center <= h - 2
    # The center of layer k is at 2 + 3*(k-1).
    # 2 + 3*(k-1) <= h - 2 => 3*(k-1) <= h - 4 => k-1 <= (h-4)/3 => k <= (h-4)/3 + 1
    if h < 4:
      num_layers = 0
    else:
      num_layers = math.floor((h - 4) / 3) + 1

    total_balls = 0
    for i in range(1, num_layers + 1):
        is_a_layer = (i % 2 == 1)
        if is_a_layer:
            # Layer A centers start at (2,2) with a 4cm grid
            # x_center = 2+4*i, x_center+2 <= l => 4+4i <= l => i <= (l-4)/4
            nx = math.floor((l - 4) / 4) + 1 if l >= 4 else 0
            ny = math.floor((w - 4) / 4) + 1 if w >= 4 else 0
            total_balls += nx * ny
        else:
            # Layer B centers start at (4,4) with a 4cm grid
            # x_center = 4+4i, x_center+2 <= l => 6+4i <= l => i <= (l-6)/4
            nx = math.floor((l - 6) / 4) + 1 if l >= 6 else 0
            ny = math.floor((w - 6) / 4) + 1 if w >= 6 else 0
            total_balls += nx * ny
            
    return total_balls

def find_optimal_box():
    """
    Searches for the box with the minimum surface area that can hold >= 27 balls.
    """
    initial_surface_area = 864  # for 12x12x12 box
    min_surface_area = initial_surface_area
    best_dims = (12, 12, 12)
    
    # Search range for dimensions (in cm, multiples of 0.5)
    # Range is chosen based on the idea that the new box will be somewhat cubic
    # and have dimensions around 10-13 cm.
    search_steps = [i * 0.5 for i in range(16, 29)] # from 8.0 to 14.0

    for l in search_steps:
        for w in search_steps:
            for h in search_steps:
                # To reduce redundant calculations, sort dimensions
                dims = tuple(sorted((l, w, h), reverse=True))
                
                surface_area = 2 * (dims[0] * dims[1] + dims[0] * dims[2] + dims[1] * dims[2])
                
                if surface_area >= min_surface_area:
                    continue
                
                num_balls = calculate_packed_balls(dims[0], dims[1], dims[2])
                
                if num_balls >= 27:
                    min_surface_area = surface_area
                    best_dims = dims

    # Use integers if dimensions are whole numbers for cleaner output
    l_out = int(best_dims[0]) if best_dims[0] == int(best_dims[0]) else best_dims[0]
    w_out = int(best_dims[1]) if best_dims[1] == int(best_dims[1]) else best_dims[1]
    h_out = int(best_dims[2]) if best_dims[2] == int(best_dims[2]) else best_dims[2]
    
    sa_out = int(min_surface_area) if min_surface_area == int(min_surface_area) else min_surface_area

    print(f"{sa_out}[box {l_out}x{w_out}x{h_out}]")

find_optimal_box()