import math

def solve_seismic_scanning():
    """
    Finds an optimal set of scanning locations in the Isis pyramid using a greedy algorithm.
    """
    # 1. Define constants and parameters
    BASE_SIDE = 150.0
    HEIGHT = 110.0
    BASE_HALF = BASE_SIDE / 2.0
    GRID_STEP = 0.5

    R_LONG = 20.0
    R_SHORT = 7.0

    V_PYRAMID = (1/3) * (BASE_SIDE**2) * HEIGHT
    V_LONG = (4/3) * math.pi * (R_LONG**3)
    V_SHORT = (4/3) * math.pi * (R_SHORT**3)
    
    # Pre-calculate a geometric factor for the distance check
    DIST_NORM = math.sqrt(HEIGHT**2 + BASE_HALF**2)

    # List to store the coordinates and radius of each placed scanner
    placed_spheres = []

    def is_sphere_fully_inside(cx, cy, cz, r):
        """Checks if a sphere is fully contained within the pyramid."""
        # Check distance to the base plane (z=0)
        if cz < r:
            return False
        
        # Check distance to the four slanted side planes.
        # The condition is derived from the plane equation `h*x + b_half*z - h*b_half = 0`
        # and the formula for distance from a point to a plane.
        required_margin = r * DIST_NORM
        if HEIGHT * BASE_HALF - HEIGHT * abs(cx) - BASE_HALF * cz < required_margin:
            return False
        if HEIGHT * BASE_HALF - HEIGHT * abs(cy) - BASE_HALF * cz < required_margin:
            return False
        
        return True

    def does_not_overlap(cx, cy, cz, r, spheres):
        """Checks if a new sphere overlaps with any already placed spheres."""
        for sp_x, sp_y, sp_z, sp_r in spheres:
            dist_sq = (cx - sp_x)**2 + (cy - sp_y)**2 + (cz - sp_z)**2
            radii_sum_sq = (r + sp_r)**2
            if dist_sq < radii_sum_sq:
                return False
        return True

    # 2. Greedy placement algorithm
    scanner_radii = [R_LONG, R_SHORT]
    for radius in scanner_radii:
        # Determine the valid range for the z-coordinate of the sphere's center
        cz_min_val = radius
        cz_max_val = HEIGHT - radius * (DIST_NORM / BASE_HALF)

        # Convert to grid indices
        z_start_idx = math.ceil(cz_min_val / GRID_STEP)
        z_end_idx = math.floor(cz_max_val / GRID_STEP)

        # Iterate through the grid from bottom to top
        for z_idx in range(z_start_idx, z_end_idx + 1):
            cz = z_idx * GRID_STEP
            
            # Determine the valid horizontal range for this height
            max_c_val = (HEIGHT * BASE_HALF - BASE_HALF * cz - radius * DIST_NORM) / HEIGHT
            if max_c_val < 0:
                continue
            
            c_max_idx = math.floor(max_c_val / GRID_STEP)

            # Iterate through x and y coordinates
            for x_idx in range(-c_max_idx, c_max_idx + 1):
                cx = x_idx * GRID_STEP
                for y_idx in range(-c_max_idx, c_max_idx + 1):
                    cy = y_idx * GRID_STEP
                    
                    # Check validity and non-overlap
                    if is_sphere_fully_inside(cx, cy, cz, radius):
                        if does_not_overlap(cx, cy, cz, radius, placed_spheres):
                            placed_spheres.append((cx, cy, cz, radius))
    
    # 3. Calculate and print the results
    n_long = sum(1 for s in placed_spheres if s[3] == R_LONG)
    n_short = sum(1 for s in placed_spheres if s[3] == R_SHORT)
    
    total_scanned_volume = n_long * V_LONG + n_short * V_SHORT
    coverage_ratio = (total_scanned_volume / V_PYRAMID) * 100

    print(f"{n_long}:{n_short}:{coverage_ratio:.1f}")

solve_seismic_scanning()
<<<5:48:28.7>>>