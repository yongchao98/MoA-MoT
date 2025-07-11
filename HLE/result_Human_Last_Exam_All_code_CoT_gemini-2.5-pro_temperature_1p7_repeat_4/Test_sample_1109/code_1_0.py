import math

def solve():
    """
    This script finds an optimal configuration for seismic scanners in a pyramid
    by using a greedy algorithm. It prioritizes placing large scanners first to
    minimize the total number of scans, searching from the outside-in to maximize
    packing efficiency within the constrained shape.

    The calculation is computationally intensive and may take several minutes to run.
    """

    # Pyramid and Scanner parameters
    H = 110.0
    L = 150.0
    R_long = 20.0
    R_short = 7.0
    step = 0.5

    # Pre-calculated values
    V_pyramid = (1.0/3.0) * (L**2) * H
    side_plane_const = math.sqrt((2*H)**2 + L**2)  # Constant for distance to side planes

    placed_scans = []

    def is_valid_placement(center, radius, existing_scans):
        """
        Checks if a sphere can be placed at the given center without
        violating the pyramid boundary and non-overlap constraints.
        """
        xc, yc, zc = center

        # Constraint 1: Sphere must be entirely inside the pyramid.
        # Check against bottom plane (z=0)
        if zc < radius:
            return False

        # Check against the four slanted side planes.
        # The distance from the center to each plane must be >= radius.
        # The plane equation is 2*H*x + L*z - L*H = 0.
        # As the center is inside, the term is negative, so distance is (L*H - 2*H*|x| - L*z) / const.
        required_dist_from_wall = radius * side_plane_const
        if (L*H - 2*H*max(abs(xc), abs(yc)) - L*zc) < required_dist_from_wall:
            return False

        # Constraint 2: Non-overlapping with existing scans.
        # The distance between centers must be >= sum of radii.
        for other_center, other_radius in existing_scans:
            min_dist_sq = (radius + other_radius)**2
            # Using squared distances to avoid costly square root operations
            dist_sq = sum((c1 - c2)**2 for c1, c2 in zip(center, other_center))
            if dist_sq < min_dist_sq:
                return False

        return True

    # --- Main Placement Loop ---
    # Process long-range scanners first, then short-range.
    for radius in [R_long, R_short]:
        # Determine valid Z range for this radius
        min_zc = radius
        # Calculate max zc possible (for a sphere at x=0, y=0)
        max_zc_val = H - radius * side_plane_const / L

        z_coords = [i * step for i in range(int(min_zc / step), int(max_zc_val / step) + 2)]

        for zc in z_coords:
            # Determine valid X/Y range for this zc (from outside-in)
            max_coord_val_num = L*H - L*zc - radius * side_plane_const
            if max_coord_val_num < 0:
                continue
            
            # Max possible coordinate value for x or y
            max_coord = max_coord_val_num / (2*H)
            
            # Iterate from outside-in, only for the octant where x >= y >= 0
            coord_range = [i * step for i in range(int(max_coord / step) + 1)]
            
            for yc in reversed(coord_range):
                for xc in reversed(coord_range):
                    if xc < yc: continue # Ensure we stay in the x >= y octant

                    # Generate symmetric candidate points from the seed (xc, yc, zc)
                    # Using a set handles duplicates (e.g., if xc=0, yc=0, or xc=yc)
                    symmetric_seeds = sorted(list(set([
                        (xc, yc), (-xc, yc), (xc, -yc), (-xc, -yc),
                        (yc, xc), (-yc, xc), (yc, -xc), (-yc, -xc)
                    ])))

                    newly_added_this_step = []
                    for x_seed, y_seed in symmetric_seeds:
                        candidate_center = (x_seed, y_seed, zc)

                        # Check if this candidate is valid against ALL previously placed scans
                        if is_valid_placement(candidate_center, radius, placed_scans):
                             # And check against scans added in this symmetric batch
                            is_ok_with_batch = True
                            for new_center, new_rad in newly_added_this_step:
                                dist_sq = sum((c1-c2)**2 for c1,c2 in zip(candidate_center, new_center))
                                if dist_sq < (radius + new_rad)**2:
                                    is_ok_with_batch = False
                                    break
                            if is_ok_with_batch:
                                newly_added_this_step.append((candidate_center, radius))
                    
                    if newly_added_this_step:
                        placed_scans.extend(newly_added_this_step)
                        # After adding a batch, we can break inner loops and move to next y/z,
                        # but continuing the search finds more packing opportunities.

    # --- Tally and Print Results ---
    num_long_scans = 0
    num_short_scans = 0
    total_scanned_volume = 0.0

    for _, radius in placed_scans:
        if radius == R_long:
            num_long_scans += 1
        else:
            num_short_scans += 1
        total_scanned_volume += (4.0/3.0) * math.pi * radius**3
        
    coverage_ratio = (total_scanned_volume / V_pyramid) * 100.0

    n = num_long_scans
    m = num_short_scans
    p_rounded = round(coverage_ratio, 1)

    print(f"{n}:{m}:{p_rounded}")


solve()
<<<5:144:45.3>>>