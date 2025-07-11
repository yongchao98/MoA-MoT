import math

# This script calculates an optimal placement for seismic scanners inside a square pyramid.
# It uses a greedy algorithm to maximize volume coverage while minimizing the number of scans.
# The algorithm first places the more efficient long-range scanners, then fills the remaining
# space with short-range scanners.

# --- Step 1: Define Constants ---
BASE_SIDE = 150.0
HEIGHT = 110.0
GRID_STEP = 0.5
R_LONG = 20.0
R_SHORT = 7.0

# Pre-calculated values for efficiency
PYRAMID_VOLUME = (1/3) * (BASE_SIDE**2) * HEIGHT
VOL_LONG = (4/3) * math.pi * (R_LONG**3)
VOL_SHORT = (4/3) * math.pi * (R_SHORT**3)
# The normal vector for the slanted faces has components (22, 0, 15) or (0, 22, 15).
# The magnitude is sqrt(22^2 + 15^2)
SQRT_709 = math.sqrt(22**2 + 15**2)

# --- Step 2: Containment Check ---
def is_contained(x, y, z, r):
    """Checks if a sphere of radius r centered at (x,y,z) is inside the pyramid."""
    # Condition 1: Sphere is above the base
    if z < r:
        return False
    
    # Condition 2: Distance to the four slanted faces
    # The equations of the faces are +/- 22x + 15z - 1650 = 0 and +/- 22y + 15z - 1650 = 0.
    # For a point (x,y,z) inside, the expression is negative. The distance is |expr|/sqrt(709).
    # We need this distance to be >= r. So, -expr >= r * sqrt(709).
    required_margin = r * SQRT_709
    
    # This single check covers the four slanted faces due to symmetry
    if 15 * z + 22 * max(abs(x), abs(y)) > 1650 - required_margin:
        return False
        
    return True

# --- Step 3: Overlap Check ---
def overlaps(x, y, z, r, placed_scans):
    """Checks if a new sphere overlaps with any in the placed_scans list."""
    for px, py, pz, pr in placed_scans:
        # Using squared distances to avoid costly sqrt operations
        dist_sq = (x - px)**2 + (y - py)**2 + (z - pz)**2
        min_dist = r + pr
        if dist_sq < min_dist**2:
            return True
    return False

# --- Step 4: Main Greedy Placement Algorithm ---
def find_optimal_scans():
    """
    Finds scanner placements using a greedy strategy.
    Returns a list of placed scans, where each item is (x, y, z, radius).
    """
    placed_scans = []
    
    # Process scanners in a list: long-range first, then short-range
    scanner_types = [("long", R_LONG), ("short", R_SHORT)]

    for name, r in scanner_types:
        print(f"Placing {name}-range scanners (r={r}m)...")
        # Determine the valid search range for the center's z-coordinate
        # Based on is_contained: 15*z <= 1650 - r*sqrt(709) for x=y=0
        z_max_val = (1650 - r * SQRT_709) / 15.0
        
        # Iterate through all possible grid points in a strategic order (bottom-up)
        z_min_steps = math.ceil(r / GRID_STEP)
        z_max_steps = math.floor(z_max_val / GRID_STEP)

        for iz in range(z_min_steps, z_max_steps + 1):
            z = iz * GRID_STEP
            
            # To speed up, only check for overlaps with nearby scans
            relevant_scans = [s for s in placed_scans if abs(z - s[2]) < r + s[3]]
            
            # Determine valid x,y range for this z-level
            max_xy_val = (1650 - 15*z - r * SQRT_709) / 22.0
            if max_xy_val < 0: continue
            
            xy_max_steps = math.floor(max_xy_val / GRID_STEP)
            
            # Iterate x,y from the center outwards to prioritize central placement
            ordered_xy_steps = [0] + [val for i in range(1, xy_max_steps + 1) for val in (i, -i)]
            
            for iy in ordered_xy_steps:
                y = iy * GRID_STEP
                for ix in ordered_xy_steps:
                    x = ix * GRID_STEP

                    # The loop bounds ensure containment. Now, just check for overlaps.
                    if not overlaps(x, y, z, r, relevant_scans):
                        placed_scans.append((x, y, z, r))

    return placed_scans

# --- Steps 5, 6, 7: Calculate and Print Results ---
if __name__ == "__main__":
    print("Running seismic scan placement optimization...")
    final_scans = find_optimal_scans()
    print("Calculation complete.")

    num_long = sum(1 for s in final_scans if s[3] == R_LONG)
    num_short = sum(1 for s in final_scans if s[3] == R_SHORT)
    
    total_scanned_volume = num_long * VOL_LONG + num_short * VOL_SHORT
    coverage_ratio = (total_scanned_volume / PYRAMID_VOLUME) * 100

    print("\n" + "="*40)
    print("          Optimal Scanning Results")
    print("="*40)
    print(f"Number of long-range scans (n): {num_long}")
    print(f"Number of short-range scans (m): {num_short}")
    print(f"Total number of scans: {num_long + num_short}")
    
    print("\n--- Coverage Calculation ---")
    print(f"Equation: (n * V_long + m * V_short) / V_pyramid")
    print(f"Values: ({num_long} * {VOL_LONG:.1f}m^3) + ({num_short} * {VOL_SHORT:.1f}m^3) = {total_scanned_volume:.1f} m^3")
    print(f"Coverage = {total_scanned_volume:.1f}m^3 / {PYRAMID_VOLUME:.1f}m^3 = {coverage_ratio:.1f}%")

    print("\n--- Final Answer ---")
    print("Format: n:m:p")
    print(f"{num_long}:{num_short}:{coverage_ratio:.1f}")