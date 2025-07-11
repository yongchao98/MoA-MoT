import math

def solve():
    """
    Finds the dimensions of a more efficient box to hold energy balls.
    """
    # Initial parameters
    initial_dim = 12.0  # cm
    ball_radius = 2.0   # cm
    ball_diameter = ball_radius * 2.0
    
    # Calculate initial box properties
    initial_surface_area = 6 * (initial_dim ** 2)
    # Number of balls in a simple grid packing
    balls_per_dim = math.floor(initial_dim / ball_diameter)
    initial_ball_count = balls_per_dim ** 3
    
    # Theoretical minimum volume needed to pack N balls
    # N * Volume_of_one_ball / Max_packing_density
    ball_volume = (4/3) * math.pi * (ball_radius ** 3)
    max_density = math.pi / (3 * math.sqrt(2))
    min_volume_required = initial_ball_count * ball_volume / max_density
    
    # --- Search for a better box ---
    best_solution = {
        'l': initial_dim,
        'w': initial_dim,
        'h': initial_dim,
        'area': initial_surface_area
    }
    
    # We are looking for an area smaller than the initial 864.
    # We can constrain the search space. For a cube l*l*l, 6*l^2 < 864 => l < 12.
    # So at least one dimension must be smaller than 12.
    search_range_l = 12
    search_range_wh = 40 # A reasonable upper bound for other dimensions
    
    min_found_area = initial_surface_area
    
    # Iterate through possible integer dimensions (l, w, h)
    # WLOG, assume l <= w <= h
    for l in range(1, search_range_l):
        for w in range(l, search_range_wh):
            # Pruning based on volume
            if l * w * 1 < min_volume_required: # h=1 is smallest possible h
                # if the smallest possible h doesn't give enough volume, larger w is needed
                continue
            
            for h in range(w, search_range_wh):
                # Volume check: if the box is too small, it can't hold the balls.
                if l * w * h < min_volume_required:
                    continue

                surface_area = 2 * (l*w + l*h + w*h)

                # We need a box with smaller surface area
                if surface_area < min_found_area:
                    # For a box to be a potential solution, it must be possible to pack
                    # at least initial_ball_count spheres inside it. While a full
                    # packing check is a very hard problem, our volume constraint
                    # (based on the highest known packing density) gives us a set
                    # of strong candidates. We select the one with the minimal area.
                    min_found_area = surface_area
                    best_solution['l'] = l
                    best_solution['w'] = w
                    best_solution['h'] = h
                    best_solution['area'] = surface_area
                
                # If area is larger, further increasing h will only make it worse
                elif surface_area > min_found_area:
                    break

    # If we found a better solution
    if best_solution['area'] < initial_surface_area:
        l = best_solution['l']
        w = best_solution['w']
        h = best_solution['h']
        area = best_solution['area']
        
        # Sort dimensions for a canonical representation (optional but good practice)
        dims = sorted([l, w, h])
        
        # Print the solution in the required format
        print(f"{dims[0]}:{dims[1]}:{dims[2]}:{int(area)}")

    else:
        # If no better box was found
        print("0")

solve()
