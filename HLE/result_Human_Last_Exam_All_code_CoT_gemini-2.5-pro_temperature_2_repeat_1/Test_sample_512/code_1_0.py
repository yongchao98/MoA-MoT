import math

def solve():
    """
    This function searches for a more material-efficient container for 27 energy balls.
    """
    # --- Initial Cube Details ---
    initial_side = 12.0
    initial_surface_area = 6 * initial_side**2
    
    ball_radius = 2.0
    ball_diameter = 4.0
    num_balls = 27

    # --- Plausibility Check based on Packing Density ---
    # Maximum known packing density (Kepler conjecture) is for FCC or HCP lattices.
    max_packing_density = math.pi / (3 * math.sqrt(2))  # ~0.74048
    volume_of_all_balls = num_balls * (4/3) * math.pi * ball_radius**3
    
    # The container must have enough volume to hold the balls, even with the densest packing.
    min_plausible_volume = volume_of_all_balls / max_packing_density

    # --- Search for the Optimal Cuboid Container ---
    min_sa = initial_surface_area
    best_box_dims = None
    
    # Set a search range for box dimensions (in cm). 
    # The optimal box should have dimensions "around" the original 12 cm.
    search_min_dim = 4.0 # A side must be at least the ball diameter.
    search_max_dim = 15.0 # Unlikely to be much larger than the original side.
    step = 0.5 # Precision constraint.

    # Iterate through possible dimensions for L, W, H
    l_dim = search_min_dim
    while l_dim <= search_max_dim:
        w_dim = l_dim # Assume L <= W <= H to avoid duplicate checks
        while w_dim <= search_max_dim:
            h_dim = w_dim
            while h_dim <= search_max_dim:
                
                volume = l_dim * w_dim * h_dim
                # Check if the volume is sufficient
                if volume < min_plausible_volume:
                    h_dim += step
                    continue

                surface_area = 2 * (l_dim * w_dim + l_dim * h_dim + w_dim * h_dim)

                if surface_area < min_sa:
                    min_sa = surface_area
                    best_box_dims = (l_dim, w_dim, h_dim)
                
                h_dim += step
            w_dim += step
        l_dim += step

    if best_box_dims:
        l, w, h = best_box_dims
        # Format for clean integer output if .0
        l_str = str(int(l)) if l == int(l) else str(l)
        w_str = str(int(w)) if w == int(w) else str(w)
        h_str = str(int(h)) if h == int(h) else str(h)
        sa_str = str(int(min_sa)) if min_sa == int(min_sa) else str(min_sa)
        
        print(f"{sa_str}[box {l_str}x{w_str}x{h_str}]")
    else:
        print(0)

solve()