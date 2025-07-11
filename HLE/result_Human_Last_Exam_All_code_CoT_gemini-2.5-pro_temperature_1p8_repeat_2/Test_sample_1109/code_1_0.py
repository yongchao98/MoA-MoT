import math
import sys

def solve_scanning_problem():
    """
    This script finds an optimal placement for seismic scanners inside a pyramid
    using a greedy algorithm. It first places larger scanners, then fills gaps
    with smaller ones, prioritizing locations high up and near the center.
    The final output is the count of each scanner type and the total volume coverage.
    """
    
    # Define pyramid and scanner constants
    BASE_SIDE = 150.0
    HEIGHT = 110.0
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5
    SQRT_709 = math.sqrt(22**2 + 15**2)

    # --- Helper Functions ---

    def is_contained(center, radius, constant_1650):
        """Checks if a sphere is fully contained within the pyramid."""
        cx, cy, cz = center
        
        # Check against the base plane (z=0)
        if cz < radius:
            return False
            
        # Check against the four slanted side planes
        required_dist_scaled = radius * SQRT_709
        
        # Check X-axis constraint
        if 22 * abs(cx) + 15 * cz + required_dist_scaled > constant_1650:
            return False
        
        # Check Y-axis constraint
        if 22 * abs(cy) + 15 * cz + required_dist_scaled > constant_1650:
            return False
            
        return True

    def is_overlapping(center, radius, placed_scanners):
        """Checks if a new scanner overlaps with any already placed scanners."""
        for scanner in placed_scanners:
            p_center = scanner['center']
            p_radius = scanner['radius']
            
            # Using squared distances is faster as it avoids the square root
            dist_sq = (center[0] - p_center[0])**2 + \
                      (center[1] - p_center[1])**2 + \
                      (center[2] - p_center[2])**2
                      
            min_dist_sq = (radius + p_radius)**2
            if dist_sq < min_dist_sq:
                return True
                
        return False

    # --- Main Logic ---

    print("Starting the search for optimal scanner placement...")
    print("This process simulates placing scanners one by one and may take a few moments.")

    placed_scanners = []
    
    # Derived from pyramid geometry (22 * 75 or 15 * 110)
    # Using a small epsilon to handle floating point inaccuracies
    pyramid_constant = 1650.0001 
    
    # --- Phase 1: Place Long-Range Scanners ---
    
    print("\nPhase 1: Placing long-range (20m) scanners...")
    
    # Define search space for long-range scanners
    max_cz_long = int((pyramid_constant / 15) / GRID_STEP) * GRID_STEP
    min_cz_long = R_LONG
    
    z_steps_long = int((max_cz_long - min_cz_long) / GRID_STEP) + 1
    
    # Sort candidates to prioritize high-z and central locations
    candidate_centers = []
    for i in range(z_steps_long):
        cz = max_cz_long - i * GRID_STEP
        # Determine max |x| and |y| for this z-level to create a bounding box
        max_c = (pyramid_constant - 15 * cz - R_LONG * SQRT_709) / 22
        if max_c < 0: continue
        max_c_steps = int(max_c / GRID_STEP)
        
        for j in range(-max_c_steps, max_c_steps + 1):
            for k in range(-max_c_steps, max_c_steps + 1):
                # We only need to check the candidate's containment once
                center_candidate = (k * GRID_STEP, j * GRID_STEP, cz)
                if is_contained(center_candidate, R_LONG, pyramid_constant):
                    candidate_centers.append(center_candidate)

    # Sort candidates: high z -> close to yz-plane -> close to xz-plane
    candidate_centers.sort(key=lambda p: (-p[2], abs(p[0]), abs(p[1])))

    for center in candidate_centers:
        if not is_overlapping(center, R_LONG, placed_scanners):
            placed_scanners.append({'center': center, 'radius': R_LONG})

    n_long = len(placed_scanners)
    print(f"Placed {n_long} long-range scanners.")
    
    # --- Phase 2: Place Short-Range Scanners ---

    print("\nPhase 2: Filling gaps with short-range (7m) scanners...")

    # Define search space for short-range scanners
    max_cz_short = int((pyramid_constant / 15) / GRID_STEP) * GRID_STEP
    min_cz_short = R_SHORT
    z_steps_short = int((max_cz_short - min_cz_short) / GRID_STEP) + 1
    
    # Generate and sort candidates
    candidate_centers_short = []
    for i in range(z_steps_short):
        cz = max_cz_short - i * GRID_STEP
        max_c = (pyramid_constant - 15 * cz - R_SHORT * SQRT_709) / 22
        if max_c < 0: continue
        max_c_steps = int(max_c / GRID_STEP)
        
        for j in range(-max_c_steps, max_c_steps + 1):
             for k in range(-max_c_steps, max_c_steps + 1):
                center_candidate = (k * GRID_STEP, j * GRID_STEP, cz)
                if is_contained(center_candidate, R_SHORT, pyramid_constant):
                    candidate_centers_short.append(center_candidate)

    candidate_centers_short.sort(key=lambda p: (-p[2], abs(p[0]), abs(p[1])))

    for center in candidate_centers_short:
        if not is_overlapping(center, R_SHORT, placed_scanners):
             placed_scanners.append({'center': center, 'radius': R_SHORT})

    n_short = len(placed_scanners) - n_long
    print(f"Placed {n_short} short-range scanners.")
    
    # --- Phase 3: Calculate and Print Results ---
    PYRAMID_VOLUME = (1/3) * BASE_SIDE**2 * HEIGHT
    V_LONG = (4/3) * math.pi * R_LONG**3
    V_SHORT = (4/3) * math.pi * R_SHORT**3
    
    total_scanned_volume = (n_long * V_LONG) + (n_short * V_SHORT)
    coverage_ratio = (total_scanned_volume / PYRAMID_VOLUME) * 100
    
    print("\n--- Calculation Breakdown ---")
    print(f"Total Pyramid Volume: {PYRAMID_VOLUME:.1f} m^3")
    print(f"Number of long-range scans (n): {n_long}")
    print(f"Number of short-range scans (m): {n_short}")
    print(f"Total scanned volume: ({n_long} * {V_LONG:.1f}) + ({n_short} * {V_SHORT:.1f}) = {total_scanned_volume:.1f} m^3")
    print(f"Coverage ratio (p): {total_scanned_volume:.1f} / {PYRAMID_VOLUME:.1f} = {coverage_ratio:.1f}%")
    
    final_answer = f"{n_long}:{n_short}:{coverage_ratio:.1f}"
    print("\nFinal Answer (n:m:p):")
    print(final_answer)
    # The platform will capture the answer in the special brackets below
    print(f"\n<<<{final_answer}>>>")

if __name__ == '__main__':
    solve_scanning_problem()