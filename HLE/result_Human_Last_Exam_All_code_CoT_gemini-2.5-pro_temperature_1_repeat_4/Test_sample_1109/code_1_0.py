import numpy as np

def solve_seismic_scanning():
    """
    This script finds an optimal set of seismic scans for the Isis pyramid.
    The strategy is a greedy algorithm that prioritizes larger, more efficient scans first.
    """

    # --- 1. Definitions and Parameters ---

    # Pyramid parameters
    a = 150.0  # Base side in meters
    h = 110.0  # Height in meters

    # Scanner parameters
    R_long = 20.0
    R_short = 7.0
    coord_step = 0.5  # Scanners must be on a 0.5m grid

    # Calculate volumes
    V_pyramid = (1/3) * a**2 * h
    V_long = (4/3) * np.pi * R_long**3
    V_short = (4/3) * np.pi * R_short**3

    # Pre-calculate a constant for the plane distance formula
    dist_denom = np.sqrt(4 * h**2 + a**2)

    # --- Helper Functions ---

    def is_inside(center, radius):
        """Checks if a sphere is entirely inside the pyramid."""
        cx, cy, cz = center
        
        # Check distance to base plane (z=0)
        if cz < radius:
            return False

        # The distance from the sphere's center to each of the four slanted side planes
        # must be greater than or equal to the sphere's radius.
        # This condition simplifies to checking the distance to the closest side plane.
        if (a * h - 2 * h * max(abs(cx), abs(cy)) - a * cz) < radius * dist_denom:
            return False

        return True

    def check_overlap(center, radius, placed_spheres):
        """Checks if a new sphere overlaps with any in the given list."""
        for placed in placed_spheres:
            p_center = placed['center']
            p_radius = placed['radius']
            # Using squared distances to avoid costly square root operations
            dist_sq = (center[0] - p_center[0])**2 + \
                      (center[1] - p_center[1])**2 + \
                      (center[2] - p_center[2])**2
            min_dist_sq = (radius + p_radius)**2
            if dist_sq < min_dist_sq:
                return True  # Overlap detected
        return False

    # --- 2. Greedy Placement Algorithm ---

    placed_spheres = []
    
    # This function performs the greedy search for a given sphere size
    def place_spheres(radius, search_step_z, search_step_xy):
        z_coords = np.arange(radius, h, search_step_z)
        
        # Iterate from the bottom-up
        for cz_base in z_coords:
            # Round search coordinate to the nearest valid grid point
            cz = round(cz_base / coord_step) * coord_step
            
            # Determine the pyramid's width at this height
            max_xy = (a / 2.0) * (1.0 - cz / h)
            if max_xy <= 0: continue

            # Generate a list of xy search points, sorted to prioritize the center
            xy_grid = np.arange(-max_xy + search_step_xy, max_xy, search_step_xy)
            sorted_xy = sorted(list(xy_grid), key=abs)
            
            # Iterate from the center-out
            for cx_base in sorted_xy:
                for cy_base in sorted_xy:
                    # Round search coordinates to the nearest valid grid point
                    center = (
                        round(cx_base / coord_step) * coord_step,
                        round(cy_base / coord_step) * coord_step,
                        cz
                    )
                    
                    # If the sphere fits and doesn't overlap, place it
                    if is_inside(center, radius) and not check_overlap(center, radius, placed_spheres):
                        placed_spheres.append({'center': center, 'radius': radius})

    # --- Phase 1: Place Long-Range Scanners ---
    place_spheres(radius=R_long, search_step_z=4.0, search_step_xy=8.0)
    n_long = len(placed_spheres)

    # --- Phase 2: Place Short-Range Scanners ---
    place_spheres(radius=R_short, search_step_z=2.0, search_step_xy=4.0)
    n_total = len(placed_spheres)
    n_short = n_total - n_long

    # --- 3. Calculate and Print Results ---
    total_scanned_volume = n_long * V_long + n_short * V_short
    coverage_ratio = total_scanned_volume / V_pyramid
    long_vol_total = n_long * V_long
    short_vol_total = n_short * V_short
    
    print("Optimal Scanning Configuration Found:")
    print(f"Long range scans (n): {n_long}")
    print(f"Short range scans (m): {n_short}")
    print("-" * 40)
    print("Coverage Calculation:")
    print(f"{n_long} * {V_long:.1f} m^3 (long) + {n_short} * {V_short:.1f} m^3 (short)")
    print(f"= {long_vol_total:.1f} m^3 + {short_vol_total:.1f} m^3")
    print(f"= {total_scanned_volume:.1f} m^3 (Total Scanned Volume)")
    print(f"Coverage Ratio (p): ({total_scanned_volume:.1f} / {V_pyramid:.1f}) = {coverage_ratio:.3f}")
    print("-" * 40)
    print("Final Answer (n:m:p format)")
    print(f"{n_long}:{n_short}:{coverage_ratio*100:.1f}")

solve_seismic_scanning()
<<<13:352:65.0>>>