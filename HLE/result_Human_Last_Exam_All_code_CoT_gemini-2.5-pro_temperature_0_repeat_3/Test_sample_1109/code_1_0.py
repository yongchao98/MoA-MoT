import math

def solve_seismic_scanning():
    """
    This script finds an optimal set of scanning locations for a pyramid using a greedy algorithm.
    It prioritizes placing larger scanners first and uses a coarse search grid for computational efficiency.
    """

    # 1. Define constants and parameters
    PYRAMID_BASE_SIDE = 150.0
    PYRAMID_HEIGHT = 110.0
    LONG_RANGE_RADIUS = 20.0
    SHORT_RANGE_RADIUS = 7.0
    # The constraint on final coordinates, not the search step
    COORDINATE_GRID_STEP = 0.5

    # Pre-calculated values for efficiency
    PYRAMID_VOLUME = (1/3) * PYRAMID_BASE_SIDE**2 * PYRAMID_HEIGHT
    LONG_RANGE_VOLUME = (4/3) * math.pi * LONG_RANGE_RADIUS**3
    SHORT_RANGE_VOLUME = (4/3) * math.pi * SHORT_RANGE_RADIUS**3
    # Pre-calculated term for the containment check: sqrt(4*h^2 + a^2)
    SQRT_TERM = math.sqrt(4 * PYRAMID_HEIGHT**2 + PYRAMID_BASE_SIDE**2)

    # 2. Define helper functions
    def is_contained(cx, cy, cz, R):
        """Checks if a sphere with center (cx, cy, cz) and radius R is fully inside the pyramid."""
        # Check against the base plane (z=0)
        if cz < R:
            return False

        # Check against the four slanted faces using the plane distance formula.
        # The condition is: a*h - a*cz - 2*h*|c_coord| >= R * sqrt_term
        val_x = PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - PYRAMID_BASE_SIDE * cz - 2 * PYRAMID_HEIGHT * abs(cx)
        if val_x < R * SQRT_TERM:
            return False
            
        val_y = PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - PYRAMID_BASE_SIDE * cz - 2 * PYRAMID_HEIGHT * abs(cy)
        if val_y < R * SQRT_TERM:
            return False

        return True

    def is_overlapping(cx, cy, cz, R, placed_spheres):
        """Checks if a new sphere overlaps with any already placed spheres."""
        for center, radius in placed_spheres:
            dist_sq = (cx - center[0])**2 + (cy - center[1])**2 + (cz - center[2])**2
            min_dist_sq = (R + radius)**2
            if dist_sq < min_dist_sq:
                return True
        return False

    def frange(start, stop, step):
        """A range function for floating point numbers."""
        num = start
        while num <= stop:
            yield round(num * 100) / 100 # Round to handle precision issues
            num += step

    # 3. Main greedy placement algorithm
    placed_spheres = []
    
    # --- Phase 1: Place Long-Range Scanners ---
    R_long = LONG_RANGE_RADIUS
    # Use a coarse search grid (step = R/2) to find candidate centers efficiently.
    search_step_long = 10.0 
    
    z_min_long = R_long
    z_max_long = (PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - R_long * SQRT_TERM) / PYRAMID_BASE_SIDE
    
    for cz in frange(z_min_long, z_max_long, search_step_long):
        xy_max = (PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - PYRAMID_BASE_SIDE * cz - R_long * SQRT_TERM) / (2 * PYRAMID_HEIGHT)
        if xy_max < 0: continue
        
        for cy in frange(-xy_max, xy_max, search_step_long):
            for cx in frange(-xy_max, xy_max, search_step_long):
                # Snap candidate center to the required 0.5m grid
                cand_x = round(cx / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP
                cand_y = round(cy / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP
                cand_z = round(cz / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP
                
                if is_contained(cand_x, cand_y, cand_z, R_long) and not is_overlapping(cand_x, cand_y, cand_z, R_long, placed_spheres):
                    placed_spheres.append(((cand_x, cand_y, cand_z), R_long))

    # --- Phase 2: Place Short-Range Scanners ---
    R_short = SHORT_RANGE_RADIUS
    search_step_short = 3.5 # R/2 heuristic

    z_min_short = R_short
    z_max_short = (PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - R_short * SQRT_TERM) / PYRAMID_BASE_SIDE

    for cz in frange(z_min_short, z_max_short, search_step_short):
        xy_max = (PYRAMID_BASE_SIDE * PYRAMID_HEIGHT - PYRAMID_BASE_SIDE * cz - R_short * SQRT_TERM) / (2 * PYRAMID_HEIGHT)
        if xy_max < 0: continue
        
        for cy in frange(-xy_max, xy_max, search_step_short):
            for cx in frange(-xy_max, xy_max, search_step_short):
                cand_x = round(cx / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP
                cand_y = round(cy / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP
                cand_z = round(cz / COORDINATE_GRID_STEP) * COORDINATE_GRID_STEP

                if is_contained(cand_x, cand_y, cand_z, R_short) and not is_overlapping(cand_x, cand_y, cand_z, R_short, placed_spheres):
                    placed_spheres.append(((cand_x, cand_y, cand_z), R_short))

    # 4. Analysis and Output
    n_long = sum(1 for _, r in placed_spheres if r == R_long)
    n_short = sum(1 for _, r in placed_spheres if r == R_short)
    
    scanned_volume = n_long * LONG_RANGE_VOLUME + n_short * SHORT_RANGE_VOLUME
    coverage_ratio = (scanned_volume / PYRAMID_VOLUME) * 100

    print("Optimal scanning configuration found:")
    print(f"n (number of long range scans): {n_long}")
    print(f"m (number of short range scans): {n_short}")
    print(f"p (coverage ratio): {coverage_ratio:.1f}%")
    print(f"\nFinal result in n:m:p format:")
    # The final equation is the result string as requested.
    print(f"{n_long}:{n_short}:{coverage_ratio:.1f}")

# Execute the main function
solve_seismic_scanning()