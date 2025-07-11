import math

def solve_seismic_scanning():
    """
    This script finds an optimal placement of seismic scanners in the Isis pyramid
    and calculates the resulting coverage.
    """

    # Pyramid and Scanner constants
    L = 150.0  # Pyramid base side in meters
    H = 110.0  # Pyramid height in meters
    R_LONG = 20.0
    R_SHORT = 7.0
    COORD_STEP = 0.5

    # Pre-calculated values for efficiency
    V_PYRAMID = (1/3) * L**2 * H
    V_LONG = (4/3) * math.pi * R_LONG**3
    V_SHORT = (4/3) * math.pi * R_SHORT**3
    
    # Pre-calculated normal vector magnitude for pyramid side planes
    # The condition for a sphere (center c, radius R) to be inside the pyramid is
    # derived from the distance to the pyramid's planes. For side planes, this simplifies to:
    # 110 * |cx| + 75 * cz <= 8250 - R * NORM_MAG
    NORM_MAG = math.sqrt(110**2 + 75**2)

    placed_scanners = []

    def is_valid(center, radius):
        """
        Checks if a new scanner placement is valid.
        A placement is valid if the sphere is inside the pyramid and does not overlap
        with any previously placed spheres.
        """
        cx, cy, cz = center

        # 1. Check if the sphere is fully inside the pyramid's boundaries
        # Check against the base plane (z=0)
        if cz < radius:
            return False
        
        # Check against the four side planes
        limit = 8250 - radius * NORM_MAG
        if 110 * abs(cx) + 75 * cz > limit or 110 * abs(cy) + 75 * cz > limit:
            return False

        # 2. Check for overlap with all existing scanners
        for scanner in placed_scanners:
            ex_center = scanner['center']
            ex_radius = scanner['radius']
            # Using squared distances to avoid costly square root operations
            dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(center, ex_center))
            min_dist_sq = (radius + ex_radius)**2
            if dist_sq < min_dist_sq - 1e-9:  # Use a small tolerance for floating point precision
                return False
                
        return True

    # --- Phase 1: Place Long-Range Scanners ---
    n_long = 0
    R = R_LONG
    # A search grid step that balances computation time and packing density.
    grid_search_step = 5.0

    cz_min = R
    cz_max = (8250 - R * NORM_MAG) / 75
    
    # Iterate from bottom-up (z-axis)
    z_coords = [i * grid_search_step for i in range(math.ceil(cz_min / grid_search_step), int(cz_max / grid_search_step) + 1)]

    for cz_raw in z_coords:
        cxy_max = (8250 - R * NORM_MAG - 75 * cz_raw) / 110
        if cxy_max < 0: continue
        
        # Iterate from center-out (xy-plane)
        xy_coords = [i * grid_search_step for i in range(int(cxy_max / grid_search_step) + 1)]
        xy_coords_sorted = sorted(list(set(xy_coords + [-c for c in xy_coords])), key=abs)

        for cy_raw in xy_coords_sorted:
            for cx_raw in xy_coords_sorted:
                # Snap candidate center to the required 0.5m coordinate grid
                center = (
                    round(cx_raw / COORD_STEP) * COORD_STEP,
                    round(cy_raw / COORD_STEP) * COORD_STEP,
                    round(cz_raw / COORD_STEP) * COORD_STEP
                )

                if is_valid(center, R):
                    placed_scanners.append({'center': center, 'radius': R})
                    n_long += 1

    # --- Phase 2: Place Short-Range Scanners ---
    m_short = 0
    R = R_SHORT
    # Use a finer grid for the smaller, more numerous spheres.
    grid_search_step = 2.0

    cz_min = R
    cz_max = (8250 - R * NORM_MAG) / 75
    
    z_coords = [i * grid_search_step for i in range(math.ceil(cz_min / grid_search_step), int(cz_max / grid_search_step) + 1)]

    for cz_raw in z_coords:
        cxy_max = (8250 - R * NORM_MAG - 75 * cz_raw) / 110
        if cxy_max < 0: continue
        
        xy_coords = [i * grid_search_step for i in range(int(cxy_max / grid_search_step) + 1)]
        xy_coords_sorted = sorted(list(set(xy_coords + [-c for c in xy_coords])), key=abs)

        for cy_raw in xy_coords_sorted:
            for cx_raw in xy_coords_sorted:
                center = (
                    round(cx_raw / COORD_STEP) * COORD_STEP,
                    round(cy_raw / COORD_STEP) * COORD_STEP,
                    round(cz_raw / COORD_STEP) * COORD_STEP
                )

                if is_valid(center, R):
                    placed_scanners.append({'center': center, 'radius': R})
                    m_short += 1
                    
    # --- Calculate and Print Final Result ---
    total_scanned_volume = n_long * V_LONG + m_short * V_SHORT
    coverage_ratio = (total_scanned_volume / V_PYRAMID) * 100
    
    # Output the number of long-range scans, short-range scans, and the coverage percentage
    n_str = str(n_long)
    m_str = str(m_short)
    p_str = f"{coverage_ratio:.1f}"

    print(f"{n_str}:{m_str}:{p_str}")

solve_seismic_scanning()
<<<2:238:43.7>>>