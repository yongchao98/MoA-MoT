import math

def solve_seismic_scanning():
    """
    This script implements a greedy algorithm to find an optimal set of scanning locations
    inside a square pyramid. It prioritizes placing larger, long-range scans first
    to maximize volume coverage efficiently, then fills remaining gaps with short-range scans.
    The search for locations is performed systematically on a 0.5m grid, from the bottom-up
    and center-out, ensuring a deterministic and well-packed result.
    """

    # 1. Define constants and pyramid geometry
    a = 150.0  # Base side in meters
    h = 110.0  # Height in meters
    R_long = 20.0
    R_short = 7.0
    step = 0.5

    # Pre-calculate geometric properties of the pyramid
    V_pyramid = (1.0 / 3.0) * a**2 * h
    k = 2.0 * h / a  # Proportionality constant for pyramid faces
    k_inv = a / (2.0 * h)
    sqrt_term = math.sqrt(1.0 + k**2)

    def get_safe_zone_params(R):
        """Calculates the parameters of the safe inner pyramid for sphere centers."""
        C_R = R * sqrt_term
        z_min = R
        z_max = h - C_R
        return z_min, z_max, C_R

    def is_overlapping(center, R, spheres):
        """Checks if a new sphere overlaps with a list of existing spheres."""
        cx, cy, cz = center
        # Use squared distances to avoid costly sqrt operations and handle floating point issues.
        # A small tolerance (1e-9) is subtracted to prevent issues with spheres that are just touching.
        min_dist_sq_threshold = (R + R)**2 - 1e-9
        
        for s_center, s_R in spheres:
            px, py, pz = s_center
            dist_sq = (cx - px)**2 + (cy - py)**2 + (cz - pz)**2
            if dist_sq < (R + s_R)**2 - 1e-9:
                return True
        return False

    def pack_spheres(R, existing_spheres):
        """
        Greedily places spheres of a given radius within the pyramid's safe zone.
        """
        newly_placed = []
        # The list of all spheres to check against includes existing ones and newly placed ones.
        all_spheres = existing_spheres + newly_placed
        
        z_min, z_max, C_R = get_safe_zone_params(R)
        
        # Iterate from the bottom of the safe zone upwards.
        z_coords = [z_min + i * step for i in range(int((z_max - z_min) / step) + 1)]

        for cz in z_coords:
            # Determine the maximum valid x/y coordinate for the center at this height.
            max_cxy = (h - cz - C_R) * k_inv
            if max_cxy < 0:
                continue

            # Iterate from the center outwards for x and y coordinates.
            xy_bound = int(max_cxy / step)
            xy_coords = [i * step for i in range(-xy_bound, xy_bound + 1)]
            sorted_xy = sorted(xy_coords, key=abs)

            for cy in sorted_xy:
                if abs(cy) > max_cxy:
                    continue
                for cx in sorted_xy:
                    if abs(cx) > max_cxy:
                        continue
                    
                    center = (cx, cy, cz)
                    
                    # A sphere is placed only if it doesn't overlap with any other sphere.
                    if not is_overlapping(center, R, all_spheres):
                        newly_placed.append((center, R))
                        # Update the list of all spheres to include the one just added.
                        all_spheres.append((center, R))
                        
        return newly_placed

    # --- Main Execution ---
    
    # Place long-range spheres first
    long_range_spheres = pack_spheres(R_long, [])
    n = len(long_range_spheres)

    # Then, place short-range spheres in the remaining space
    short_range_spheres = pack_spheres(R_short, long_range_spheres)
    m = len(short_range_spheres)
    
    # Calculate final volume coverage
    V_long = n * (4.0/3.0) * math.pi * R_long**3
    V_short = m * (4.0/3.0) * math.pi * R_short**3
    V_scanned = V_long + V_short
    coverage = (V_scanned / V_pyramid) * 100

    # Print the final result in the specified format n:m:p
    print(f"{n}:{m}:{coverage:.1f}")

solve_seismic_scanning()
<<<6:112:35.3>>>