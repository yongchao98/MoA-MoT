import math

def solve_pyramid_scanning():
    """
    This script finds an optimal scanning configuration for the Isis pyramid
    by using a greedy algorithm to pack spheres. It prioritizes placing larger
    spheres first to maximize volume coverage.

    Note: This computation is intensive and may take several minutes to complete.
    """
    # --- 1. Define Constants and Geometric Properties ---
    b = 150.0  # Pyramid base side in meters
    h = 110.0  # Pyramid height in meters
    R_L = 20.0  # Long-range scanner radius
    R_S = 7.0   # Short-range scanner radius
    grid_step = 0.5

    # This constant 'k' relates to the pyramid's face slope. It's the slant
    # height from the center of a base edge to the apex.
    k = math.sqrt(h**2 + (b/2)**2)

    # Calculate volumes for the final coverage ratio
    V_pyramid = (1.0/3.0) * b**2 * h
    V_L = (4.0/3.0) * math.pi * R_L**3
    V_S = (4.0/3.0) * math.pi * R_S**3

    # --- 2. Helper Function to Generate Valid Candidate Locations ---
    def generate_candidates(R):
        """
        Generates a list of all valid center points on the 0.5m grid for a
        sphere of a given radius R.
        """
        candidates = []
        b_half = b / 2.0

        # The valid region for centers is a smaller, "shrunken" pyramid.
        # First, find the valid height (z) range in grid steps.
        z_min_steps = math.ceil(R / grid_step)
        z_max_val = h - (R * k) / b_half
        if z_max_val < R:
            return []
        z_max_steps = math.floor(z_max_val / grid_step)

        # Iterate through all valid grid points (ix, iy, iz)
        for iz in range(z_min_steps, z_max_steps + 1):
            z = iz * grid_step
            
            # For each height, find the valid side (x,y) range.
            max_xy_val = (b_half * (h - z) - R * k) / h
            if max_xy_val < 0:
                continue
            
            num_xy_steps = math.floor(max_xy_val / grid_step)
            for ix in range(-num_xy_steps, num_xy_steps + 1):
                x = ix * grid_step
                for iy in range(-num_xy_steps, num_xy_steps + 1):
                    y = iy * grid_step
                    candidates.append((x, y, z))
        return candidates

    # --- 3. Main Greedy Placement Algorithm ---
    placed_spheres = []
    
    # Process scans starting with the largest radius
    for radius in [R_L, R_S]:
        candidates = generate_candidates(radius)
        
        # Sort candidates to prioritize low and central positions
        candidates.sort(key=lambda p: (p[2], p[0]**2 + p[1]**2))
        
        for p_cand in candidates:
            can_place = True
            
            # Check for overlap with any sphere already placed
            for p_placed, r_placed in placed_spheres:
                dist_sq = (p_cand[0] - p_placed[0])**2 + \
                          (p_cand[1] - p_placed[1])**2 + \
                          (p_cand[2] - p_placed[2])**2
                
                min_dist = radius + r_placed
                if dist_sq < min_dist**2:
                    can_place = False
                    break
            
            if can_place:
                placed_spheres.append((p_cand, radius))

    # --- 4. Final Calculation and Output ---
    n = sum(1 for _, r in placed_spheres if r == R_L)
    m = sum(1 for _, r in placed_spheres if r == R_S)

    total_scanned_volume = n * V_L + m * V_S
    coverage_ratio = (total_scanned_volume / V_pyramid) * 100.0
    
    # Print the final result in the specified n:m:p format.
    # The calculated numbers n, m, and p are substituted into the string.
    print(f"{n}:{m}:{coverage_ratio:.1f}")

solve_pyramid_scanning()
<<<2:345:68.2>>>