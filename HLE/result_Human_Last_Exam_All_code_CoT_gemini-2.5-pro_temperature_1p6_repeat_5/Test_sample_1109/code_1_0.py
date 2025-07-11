import math
import numpy as np

def solve_seismic_scanning():
    """
    Finds an optimal placement of long- and short-range scanners in a pyramid
    using a greedy algorithm.
    """

    # 1. Define Geometry and Constants
    H = 110.0  # Pyramid height
    A_HALF = 75.0  # Pyramid base half-side
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5

    # Calculate volumes for final ratio
    PYRAMID_VOL = (1.0/3.0) * (A_HALF * 2)**2 * H
    SPHERE_VOL_LONG = (4.0/3.0) * math.pi * R_LONG**3
    SPHERE_VOL_SHORT = (4.0/3.0) * math.pi * R_SHORT**3

    def get_pyramid_half_width(z):
        """Calculates the pyramid's half-width at a given height z."""
        if 0 <= z <= H:
            return A_HALF * (1.0 - z / H)
        return 0

    def check_overlap(p_new, r_new, existing_scanners):
        """
        Checks if a new scanner overlaps with any existing ones.
        An optimization is included to quickly skip scanners that are too far apart on the z-axis.
        """
        for scanner in existing_scanners:
            p_old = scanner['center']
            r_old = scanner['R']
            required_dist = r_new + r_old

            # Optimization: if z-distance is already greater than the sum of radii, they can't overlap.
            if abs(p_new[2] - p_old[2]) >= required_dist:
                continue

            dist_sq = (p_new[0] - p_old[0])**2 + (p_new[1] - p_old[1])**2 + (p_new[2] - p_old[2])**2
            if dist_sq < required_dist**2:
                return True
        return False

    placed_scanners = []

    # 3. Greedy Placement Algorithm
    # This function implements the core greedy placement logic for a given scanner type.
    def place_scanners(r_scanner, existing_scanners):
        # Iterate z from bottom-up
        z_coords = np.arange(r_scanner, H, GRID_STEP)
        
        for z in z_coords:
            # Determine the valid range for scanner centers at this height
            max_center_coord = get_pyramid_half_width(z) - r_scanner
            if max_center_coord < 0:
                continue

            # Generate potential x, y coordinates, ordered from center outward
            xy_pos_coords = np.arange(0, max_center_coord + GRID_STEP, GRID_STEP)
            xy_candidates = []
            for y_val in xy_pos_coords:
                for x_val in xy_pos_coords:
                    if x_val == 0 and y_val == 0:
                        xy_candidates.append((0.0, 0.0))
                    elif y_val == 0:
                        xy_candidates.extend([(x_val, 0.0), (-x_val, 0.0)])
                    elif x_val == 0:
                         xy_candidates.extend([(0.0, y_val), (0.0, -y_val)])
                    elif x_val == y_val:
                        xy_candidates.extend([(x_val, y_val), (-x_val, y_val), (x_val, -y_val), (-x_val, -y_val)])
                    else: # x != y and both are non-zero
                        xy_candidates.extend([
                            (x_val, y_val), (-x_val, y_val), (x_val, -y_val), (-x_val, -y_val),
                            (y_val, x_val), (-y_val, x_val), (y_val, -x_val), (-y_val, -x_val)
                        ])
            
            # Remove duplicates and sort by distance from center to prioritize central placement
            unique_candidates = sorted(list(set(xy_candidates)), key=lambda p: p[0]**2 + p[1]**2)

            for x, y in unique_candidates:
                p_new = (x, y, z)
                if not check_overlap(p_new, r_scanner, existing_scanners):
                    existing_scanners.append({'center': p_new, 'R': r_scanner})
        return existing_scanners

    # --- Phase 1: Place Long-Range Scanners ---
    placed_scanners = place_scanners(R_LONG, placed_scanners)
    num_long_scans = len(placed_scanners)

    # --- Phase 2: Place Short-Range Scanners ---
    placed_scanners = place_scanners(R_SHORT, placed_scanners)
    num_short_scans = len(placed_scanners) - num_long_scans

    # 4. Calculate and Report Results
    total_scanned_vol = (num_long_scans * SPHERE_VOL_LONG) + (num_short_scans * SPHERE_VOL_SHORT)
    coverage_ratio = (total_scanned_vol / PYRAMID_VOL) * 100

    print(f"{num_long_scans}:{num_short_scans}:{coverage_ratio:.1f}%")

solve_seismic_scanning()
<<<14:152:68.3%>>>