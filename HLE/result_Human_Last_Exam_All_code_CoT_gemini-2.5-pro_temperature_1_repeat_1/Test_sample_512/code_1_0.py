import math

def solve():
    """
    Solves the container optimization problem.
    """
    ball_radius = 2.0
    ball_diameter = 2 * ball_radius
    
    # Initial container setup
    initial_side = 12.0
    initial_surface_area = 6 * initial_side**2
    initial_ball_count = (initial_side // ball_diameter)**3
    
    min_surface_area = initial_surface_area
    best_solution = "0"

    # --- Strategy 1: Search for a better Rectangular Box ---
    
    # This function calculates how many balls fit in a staggered packing
    def count_balls_staggered(L, W, H):
        if L < ball_diameter or W < ball_diameter or H < ball_diameter:
            return 0
        
        # Layer A is a simple grid
        nx_a = math.floor(L / ball_diameter)
        ny_a = math.floor(W / ball_diameter)
        balls_layer_a = nx_a * ny_a
        
        # Layer B is placed in the hollows of layer A
        # The grid for layer B is offset by ball_radius
        nx_b = math.floor((L - ball_radius) / ball_diameter)
        ny_b = math.floor((W - ball_radius) / ball_diameter)
        balls_layer_b = nx_b * ny_b

        # Vertical distance between staggered layers.
        # d_horizontal^2 = (ball_radius)^2 + (ball_radius)^2 = 8
        # d_vertical^2 = ball_diameter^2 - d_horizontal^2 = 16 - 8 = 8
        # d_vertical = sqrt(8) ~= 2.828
        # With 0.5cm precision, the minimum separation is 3.0cm
        dz = 3.0 
        
        total_balls = 0
        current_z = ball_radius
        layer_type = 'A'
        
        while True:
            # Check if the ball fits vertically
            if current_z + ball_radius > H:
                break
            
            if layer_type == 'A':
                if balls_layer_a > 0:
                    total_balls += balls_layer_a
                    current_z += dz
                    layer_type = 'B'
                else: # If layer A is empty, so is B
                    break
            else: # Layer B
                if balls_layer_b > 0:
                    total_balls += balls_layer_b
                    current_z += dz
                    layer_type = 'A'
                else: # If layer B is empty, we can't continue staggered packing
                    break
        return total_balls

    # Search range for box dimensions
    # We search around the optimal cube shape
    search_range = [x * 0.5 for x in range(int(8 * 2), int(16 * 2))]

    for l in search_range:
        for w in search_range:
            # Let's check for a promising cuboid that is almost a cube.
            # A very long box is unlikely to be efficient.
            if w > l * 1.5: continue
            for h in search_range:
                if h > w * 1.5: continue

                surface_area = 2 * (l*w + l*h + w*h)
                
                # If SA is not better, skip
                if surface_area >= min_surface_area:
                    continue
                
                # Check packing density
                # Strategy 1: Simple Grid Packing
                grid_packed_balls = (l // ball_diameter) * (w // ball_diameter) * (h // ball_diameter)

                if grid_packed_balls >= initial_ball_count:
                    min_surface_area = surface_area
                    best_solution = f"{surface_area}[box {l}x{w}x{h}]"
                    continue # Found a better grid solution, continue search for even better

                # Strategy 2: Staggered Packing
                # Try all 3 orientations for the staggered packing height
                staggered_balls_h = count_balls_staggered(l, w, h)
                staggered_balls_w = count_balls_staggered(l, h, w)
                staggered_balls_l = count_balls_staggered(w, h, l)
                
                max_staggered_balls = max(staggered_balls_h, staggered_balls_w, staggered_balls_l)

                if max_staggered_balls >= initial_ball_count:
                    min_surface_area = surface_area
                    best_solution = f"{surface_area}[box {l}x{w}x{h}]"

    # Based on extensive search, no better container is found with these packing strategies.
    # While exotic packings exist, they are unlikely to conform to the 0.5cm coordinate grid
    # without losing their efficiency advantage. The most plausible answer is that no such
    # container can be designed under the given constraints.
    # The provided search algorithm confirms this for the tested packing strategies.

    # A promising candidate that ultimately fails is a 12x12x11.5 box.
    # SA = 2*(12*12 + 12*11.5 + 12*11.5) = 840.0 < 864.
    # However, max balls it can hold with staggered packing is 22.
    # balls_in_12_12_11_5 = count_balls_staggered(12, 12, 11.5) -> 22
    
    # Another promising candidate: Cylinder with R=6.5, H=14.5
    # SA = 2*pi*R*(R+H) = 2*pi*6.5*(6.5+14.5) = 273*pi ~ 857.6 < 864
    # However, packing 7 balls per layer requires a radius of at least 7.31cm (R=7.5),
    # which makes the surface area too large.
    
    print(best_solution)

solve()
<<<0>>>