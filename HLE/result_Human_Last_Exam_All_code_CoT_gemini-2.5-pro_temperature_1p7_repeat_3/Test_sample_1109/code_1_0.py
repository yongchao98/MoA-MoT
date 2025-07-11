import math

def solve():
    """
    Solves the seismic scanner packing problem for the Isis pyramid.
    """
    # --- Step 1: Define Constants ---
    H = 110.0
    B_SIDE = 150.0
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5

    V_PYRAMID = (1/3) * B_SIDE**2 * H
    B_HALF = B_SIDE / 2.0
    
    # Pre-calculate a constant for the containment check to improve performance
    NORM_SLANT = math.sqrt(H**2 + B_HALF**2)
    HB_HALF_PRODUCT = H * B_HALF

    placed_scanners = []

    # --- Step 2: Define Helper Function for Validation ---
    def is_valid(cx, cy, cz, r, existing_scanners):
        """
        Checks if a new scanner placement is valid (within pyramid and not overlapping).
        """
        # Constraint 1: Containment
        # Check if the sphere is above the base plane by at least its radius
        if cz < r:
            return False
        # Check if the sphere is inside the four slant faces
        if H * max(abs(cx), abs(cy)) + B_HALF * cz + r * NORM_SLANT > HB_HALF_PRODUCT:
            return False

        # Constraint 2: Non-overlapping with existing scanners
        for px, py, pz, pr in existing_scanners:
            dist_sq = (cx - px)**2 + (cy - py)**2 + (cz - pz)**2
            if dist_sq < (r + pr)**2:
                return False
        return True

    # --- Step 3: Greedy Packing Algorithm ---
    def pack_scanners(radius, scanners_list):
        """
        Iterates through possible locations to place scanners of a given radius.
        """
        r = radius
        # Determine the iteration range for the center's z-coordinate
        cz_min = r
        # The highest possible z-coordinate for a scanner's center
        cz_max = (HB_HALF_PRODUCT - r * NORM_SLANT) / B_HALF

        z = cz_min
        while z <= cz_max:
            # For this z-layer, find the maximum possible x or y coordinate for a center
            xy_limit = (HB_HALF_PRODUCT - B_HALF * z - r * NORM_SLANT) / H
            if xy_limit < 0:
                z += GRID_STEP
                continue

            # Generate a list of candidate (cx, cy) points for this layer
            layer_candidates = []
            max_coord = math.floor(xy_limit / GRID_STEP) * GRID_STEP
            x = -max_coord
            while x <= max_coord:
                y = -max_coord
                while y <= max_coord:
                    # A quick circular bound check to prune some candidates
                    if math.sqrt(x**2 + y**2) <= xy_limit + GRID_STEP:
                       layer_candidates.append((x, y))
                    y += GRID_STEP
                x += GRID_STEP
            
            # Sort candidates from the outside in (descending distance from Z-axis)
            # This heuristic helps achieve a denser packing.
            layer_candidates.sort(key=lambda p: p[0]**2 + p[1]**2, reverse=True)

            for cx, cy in layer_candidates:
                if is_valid(cx, cy, z, r, scanners_list):
                    scanners_list.append((cx, cy, z, r))
            
            z += GRID_STEP

    # --- Step 4: Execute Packing and Calculate Results ---
    
    # Phase 1: Pack long-range scanners
    pack_scanners(R_LONG, placed_scanners)
    n = len(placed_scanners)

    # Phase 2: Pack short-range scanners in the remaining gaps
    pack_scanners(R_SHORT, placed_scanners)
    m = len(placed_scanners) - n

    # Calculate total volume and coverage ratio
    v_long_sphere = (4/3) * math.pi * R_LONG**3
    v_short_sphere = (4/3) * math.pi * R_SHORT**3
    
    total_volume_scanned = n * v_long_sphere + m * v_short_sphere
    coverage_ratio = (total_volume_scanned / V_PYRAMID) * 100

    print(f"{n}:{m}:{coverage_ratio:.1f}")


solve()
<<<11:1516:34.0>>>