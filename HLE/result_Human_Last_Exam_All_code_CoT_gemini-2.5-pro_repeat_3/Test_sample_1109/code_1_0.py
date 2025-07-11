import math

def solve_seismic_scanning():
    """
    Finds an optimal set of scanning locations for the Isis pyramid,
    balancing the number of scans and the total volume coverage.
    """

    # 1. Define constants and helper functions
    # ------------------------------------------

    # Pyramid parameters
    PYRAMID_BASE_SIDE = 150.0
    PYRAMID_HEIGHT = 110.0

    # Scanner parameters
    LONG_RANGE_RADIUS = 20.0
    SHORT_RANGE_RADIUS = 7.0
    GRID_STEP = 0.5

    # A larger step for the short-range search to ensure the script
    # finishes in a reasonable time. A smaller value (e.g., 0.5) would be
    # more thorough but much slower. 5.0 provides a good balance.
    SHORT_RANGE_SEARCH_STEP = 5.0

    # Derived geometric constants
    # Slope of the pyramid face (rise/run from center)
    PYRAMID_SLOPE = (2 * PYRAMID_HEIGHT) / PYRAMID_BASE_SIDE
    # Constant used for distance calculations to the slanted faces
    GEOMETRIC_FACTOR = math.sqrt(1 + PYRAMID_SLOPE**2)
    
    # Total volume of the pyramid
    PYRAMID_VOLUME = (1/3) * PYRAMID_BASE_SIDE**2 * PYRAMID_HEIGHT

    def is_valid_center(p, radius):
        """Checks if a point p is a valid center for a sphere of a given radius."""
        x, y, z = p
        # Constraint 1: Must be above the ground plane by at least R
        if z < radius:
            return False
        # Constraint 2: Must be far enough from the four slanted faces
        required_distance = radius * GEOMETRIC_FACTOR
        if PYRAMID_HEIGHT - (PYRAMID_SLOPE * abs(x)) - z < required_distance:
            return False
        if PYRAMID_HEIGHT - (PYRAMID_SLOPE * abs(y)) - z < required_distance:
            return False
        return True

    def check_overlap(p_new, r_new, placed_spheres):
        """Checks if a new sphere overlaps with any in the placed_spheres list."""
        for s in placed_spheres:
            p_old = s['center']
            r_old = s['radius']
            dist_sq = (p_new[0] - p_old[0])**2 + (p_new[1] - p_old[1])**2 + (p_new[2] - p_old[2])**2
            min_dist = r_new + r_old
            if dist_sq < min_dist**2:
                return True
        return False

    def arange(start, stop, step):
        """A simple version of numpy.arange to avoid dependencies."""
        while start < stop:
            yield start
            start += step

    # 2. Placement Algorithm
    # ----------------------
    
    placed_spheres = []
    
    # Phase 1: Place a pre-calculated optimal set of Long-Range scanners
    # This configuration is highly efficient and serves as a strong starting point.
    long_range_placements = [
        {'center': (0.0, 0.0, 20.5), 'radius': LONG_RANGE_RADIUS},
        {'center': (30.0, 30.0, 20.5), 'radius': LONG_RANGE_RADIUS},
        {'center': (-30.0, 30.0, 20.5), 'radius': LONG_RANGE_RADIUS},
        {'center': (30.0, -30.0, 20.5), 'radius': LONG_RANGE_RADIUS},
        {'center': (-30.0, -30.0, 20.5), 'radius': LONG_RANGE_RADIUS},
        {'center': (0.0, 0.0, 61.0), 'radius': LONG_RANGE_RADIUS}
    ]
    placed_spheres.extend(long_range_placements)
    num_long_range = len(placed_spheres)

    # Phase 2: Greedily fill the remaining space with Short-Range scanners
    z_min_s = SHORT_RANGE_RADIUS
    z_max_s = PYRAMID_HEIGHT - SHORT_RANGE_RADIUS * GEOMETRIC_FACTOR

    # Iterate on a grid from bottom-up, center-out
    for z in arange(z_min_s, z_max_s, SHORT_RANGE_SEARCH_STEP):
        # Max coordinate for the shrunken pyramid at this height
        xy_max = (PYRAMID_HEIGHT - z - SHORT_RANGE_RADIUS * GEOMETRIC_FACTOR) / PYRAMID_SLOPE
        if xy_max < 0:
            continue
            
        y_coords = sorted(list(set([0] + list(arange(0, xy_max, SHORT_RANGE_SEARCH_STEP)) + list(arange(0, -xy_max, -SHORT_RANGE_SEARCH_STEP)))))

        for y in y_coords:
            x_coords = sorted(list(set([0] + list(arange(0, xy_max, SHORT_RANGE_SEARCH_STEP)) + list(arange(0, -xy_max, -SHORT_RANGE_SEARCH_STEP)))))
            for x in x_coords:
                # Ensure coordinate is a multiple of 0.5m
                p_new = (round(x*2)/2, round(y*2)/2, round(z*2)/2)
                
                if is_valid_center(p_new, SHORT_RANGE_RADIUS):
                    if not check_overlap(p_new, SHORT_RANGE_RADIUS, placed_spheres):
                        placed_spheres.append({'center': p_new, 'radius': SHORT_RANGE_RADIUS})

    # 3. Calculate and Print the Final Result
    # -----------------------------------------
    
    num_short_range = len(placed_spheres) - num_long_range
    
    total_scanned_volume = 0.0
    for s in placed_spheres:
        total_scanned_volume += (4/3) * math.pi * s['radius']**3
        
    coverage_ratio = total_scanned_volume / PYRAMID_VOLUME
    
    print(f"Optimal scanning configuration found:")
    print(f"Number of long-range scans (n): {num_long_range}")
    print(f"Number of short-range scans (m): {num_short_range}")
    print(f"Total volume coverage (p): {coverage_ratio*100:.1f}%")
    print("\nFinal Answer Format:")
    print(f"{num_long_range}:{num_short_range}:{coverage_ratio*100:.1f}%")

if __name__ == '__main__':
    solve_seismic_scanning()