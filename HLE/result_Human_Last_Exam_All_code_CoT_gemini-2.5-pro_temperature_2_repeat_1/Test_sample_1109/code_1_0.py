import numpy as np
import math

# --- Simulation Configuration ---
BASE_SIDE = 150.0  # meters
HEIGHT = 110.0   # meters
GRID_STEP = 0.5    # meters

LONG_RANGE_RADIUS = 20.0
SHORT_RANGE_RADIUS = 7.0

def run_simulation():
    """
    Main function to run the seismic scanning simulation.
    """
    print("Finding optimal scanning locations...")
    print("This process involves checking a large number of potential locations and may take a few minutes to complete.")

    placed_scans = []

    # --- Helper Functions ---

    def is_point_inside_pyramid(p):
        """Checks if a point (x, y, z) is inside the pyramid."""
        px, py, pz = p
        if not (0 <= pz <= HEIGHT):
            return False
        # Calculate the half-width of the pyramid at the point's height
        max_coord_at_z = (BASE_SIDE / 2.0) * (1.0 - pz / HEIGHT)
        if abs(px) > max_coord_at_z or abs(py) > max_coord_at_z:
            return False
        return True

    def is_sphere_valid(center, radius):
        """
        Checks if a sphere is entirely inside the pyramid by testing the 6
        cardinal points on the sphere's surface. This is a robust and fast check.
        """
        cx, cy, cz = center
        points_to_check = [
            (cx + radius, cy, cz), (cx - radius, cy, cz),
            (cx, cy + radius, cz), (cx, cy - radius, cz),
            (cx, cy, cz + radius), (cx, cy, cz - radius)
        ]
        for p in points_to_check:
            if not is_point_inside_pyramid(p):
                return False
        return True

    # --- Greedy Packing Algorithm ---

    # Process larger scanners first to maximize volume per scan
    scan_radii = [LONG_RANGE_RADIUS, SHORT_RANGE_RADIUS]

    for radius in scan_radii:
        print(f"\nPlacing scanners with radius {radius}m...")
        # The lowest a center can be is at a height equal to its own radius
        z_min = radius
        z_coords = np.arange(z_min, HEIGHT, GRID_STEP)

        # Iterate through all possible grid locations from bottom to top
        for z in z_coords:
            # Determine search bounds for this z-level to speed up computation
            max_half_width_at_z = (BASE_SIDE / 2.0) * (1.0 - z / HEIGHT)
            
            # The center of a sphere cannot be closer to the edge than its radius.
            # This provides a tighter bounding box for our search loop.
            search_limit = max_half_width_at_z - radius
            if search_limit < 0:
                continue

            limit_in_steps = int(search_limit / GRID_STEP)
            
            # Iterate through y and x coordinates (negative to positive)
            for y_step in range(-limit_in_steps, limit_in_steps + 1):
                y = y_step * GRID_STEP
                for x_step in range(-limit_in_steps, limit_in_steps + 1):
                    x = x_step * GRID_STEP
                    
                    center = (x, y, z)
                    
                    # 1. Check for overlap with already placed scans
                    is_overlapping = False
                    for placed_scan in placed_scans:
                        p_center = placed_scan['center']
                        p_radius = placed_scan['radius']
                        
                        # Use squared distances to avoid costly sqrt operations
                        dist_sq = (center[0] - p_center[0])**2 + \
                                  (center[1] - p_center[1])**2 + \
                                  (center[2] - p_center[2])**2
                        min_dist_sq = (radius + p_radius)**2
                        
                        if dist_sq < min_dist_sq:
                            is_overlapping = True
                            break
                    if is_overlapping:
                        continue
                        
                    # 2. If no overlap, check if the sphere is valid
                    if is_sphere_valid(center, radius):
                        placed_scans.append({'center': center, 'radius': radius})

    print("\n--- Calculation Results ---")
    
    # --- Analysis ---
    n = sum(1 for s in placed_scans if s['radius'] == LONG_RANGE_RADIUS)
    m = sum(1 for s in placed_scans if s['radius'] == SHORT_RANGE_RADIUS)

    vol_pyramid = (1.0/3.0) * (BASE_SIDE**2) * HEIGHT
    vol_long_scan = (4.0/3.0) * math.pi * (LONG_RANGE_RADIUS**3)
    vol_short_scan = (4.0/3.0) * math.pi * (SHORT_RANGE_RADIUS**3)
    
    vol_scanned = (n * vol_long_scan) + (m * vol_short_scan)
    coverage_ratio = (vol_scanned / vol_pyramid) * 100.0

    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print(f"Total scans: {n + m}")
    print("\nEquation for coverage ratio (p):")
    print(f"p = (n * Vol_long + m * Vol_short) / Vol_pyramid * 100")
    print(f"p = ({n} * {vol_long_scan:.2f} + {m} * {vol_short_scan:.2f}) / {vol_pyramid:.2f} * 100")
    print(f"p = ({n * vol_long_scan:.2f} + {m * vol_short_scan:.2f}) / {vol_pyramid:.2f} * 100")
    print(f"p = {vol_scanned:.2f} / {vol_pyramid:.2f} * 100")
    print(f"p = {coverage_ratio:.1f}%")

    # Final formatted answer
    final_answer_str = f"{n}:{m}:{coverage_ratio:.1f}"
    print(f"\nFinal Answer (n:m:p format): {final_answer_str}")
    
# Execute the simulation when the script is run
run_simulation()
<<<13:92:54.6>>>