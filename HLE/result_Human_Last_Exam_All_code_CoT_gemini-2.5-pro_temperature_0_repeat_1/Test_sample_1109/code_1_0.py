import math
import numpy as np

def solve_pyramid_scanning():
    """
    This function solves the seismic scanning problem by finding an optimal placement
    of long-range and short-range scanners inside the Isis pyramid.
    """

    # Step 1: Define constants for the pyramid and scanners
    a = 150.0  # Pyramid base side length
    h = 110.0  # Pyramid height
    s = a / 2.0  # Half base side length

    R_long = 20.0   # Long-range scanner radius
    R_short = 7.0   # Short-range scanner radius
    step = 0.5      # Coordinate grid step

    # Pre-calculate geometric values for efficiency
    pyramid_volume = (1/3) * (a**2) * h
    sqrt_h2_s2 = math.sqrt(h**2 + s**2)

    def get_potential_centers(radius):
        """
        Generates a list of all possible center locations for a sphere of a given
        radius that fit inside the pyramid boundaries. The list is sorted based on the
        greedy strategy: top-down (z descending), then center-out (radial distance ascending).
        """
        potential_centers = []
        
        # The highest possible z for a center is where the inscribed pyramid width is 0
        max_cz = h - radius * sqrt_h2_s2 / s
        # Generate z-coordinates from top to bottom
        z_coords = np.arange(max_cz, radius - step, -step)

        for cz in z_coords:
            # Calculate the maximum allowed |x| or |y| coordinate for a center at this height
            max_xy_coord = (s * h - radius * sqrt_h2_s2 - s * cz) / h
            if max_xy_coord < 0:
                continue
            
            centers_at_this_z = []
            # Generate coordinates in the first quadrant (x>=0, y>=0)
            xy_coords = np.arange(0, max_xy_coord + step, step)
            
            for cx in xy_coords:
                for cy in xy_coords:
                    # The loops generate a square; we only need points within the pyramid's square cross-section
                    if cx > max_xy_coord or cy > max_xy_coord:
                        continue

                    # Add points for all 4 symmetric quadrants
                    centers_at_this_z.append((cx, cy, cz))
                    if cx > 0:
                        centers_at_this_z.append((-cx, cy, cz))
                    if cy > 0:
                        centers_at_this_z.append((cx, -cy, cz))
                    if cx > 0 and cy > 0:
                        centers_at_this_z.append((-cx, -cy, cz))

            # Sort this z-level's unique centers by distance from the z-axis
            # Use a set to handle duplicates on axes, then sort
            unique_centers_at_z = sorted(list(set(centers_at_this_z)), key=lambda p: (p[0]**2 + p[1]**2))
            potential_centers.extend(unique_centers_at_z)
            
        return potential_centers

    def check_overlap(center, radius, existing_spheres):
        """Checks if a new sphere overlaps with any in a given list."""
        cx, cy, cz = center
        for (px, py, pz, pR) in existing_spheres:
            # Using squared distances to avoid costly square root operations
            dist_sq = (cx - px)**2 + (cy - py)**2 + (cz - pz)**2
            min_dist_sq = (radius + pR)**2
            if dist_sq < min_dist_sq:
                return True  # Overlaps
        return False

    # --- Main Execution Logic ---

    placed_spheres = []
    
    # Phase 1: Pack with long-range scanners
    long_range_scans = []
    potential_long = get_potential_centers(R_long)
    for center in potential_long:
        if not check_overlap(center, R_long, placed_spheres):
            sphere_info = center + (R_long,)
            placed_spheres.append(sphere_info)
            long_range_scans.append(sphere_info)

    # Phase 2: Pack with short-range scanners in the remaining space
    short_range_scans = []
    potential_short = get_potential_centers(R_short)
    for center in potential_short:
        if not check_overlap(center, R_short, placed_spheres):
            sphere_info = center + (R_short,)
            placed_spheres.append(sphere_info)
            short_range_scans.append(sphere_info)

    # Final Step: Calculate and print the results in the specified format
    n = len(long_range_scans)
    m = len(short_range_scans)

    vol_long = n * (4/3) * math.pi * (R_long**3)
    vol_short = m * (4/3) * math.pi * (R_short**3)
    total_scanned_volume = vol_long + vol_short

    coverage_ratio = (total_scanned_volume / pyramid_volume) * 100.0
    
    p_rounded = round(coverage_ratio, 1)

    # Output the final answer string
    print(f"{n}:{m}:{p_rounded}")

solve_pyramid_scanning()
<<<2:338:20.3>>>