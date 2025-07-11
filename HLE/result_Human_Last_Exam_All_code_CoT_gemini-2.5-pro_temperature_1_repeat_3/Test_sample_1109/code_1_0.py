import math

def solve():
    """
    Calculates and verifies an optimal scanning configuration for the Isis pyramid.
    """
    # Pyramid and Scanner Parameters
    BASE_SIDE = 150.0
    HEIGHT = 110.0
    R_LONG = 20.0
    R_SHORT = 7.0
    
    # --- Geometric Helper Constants ---
    # Derived from the pyramid's face plane equations (e.g., 110x + 75z - 8250 = 0)
    PLANE_DIST_FACTOR = math.sqrt(HEIGHT**2 + (BASE_SIDE / 2)**2) # sqrt(110^2 + 75^2)
    PLANE_CONSTANT = (BASE_SIDE / 2) * HEIGHT # 75 * 110 = 8250

    # --- Volume Calculations ---
    V_PYRAMID = (1/3) * BASE_SIDE**2 * HEIGHT
    V_LONG = (4/3) * math.pi * R_LONG**3
    V_SHORT = (4/3) * math.pi * R_SHORT**3

    # --- Proposed Optimal Scan Locations ---
    # n=4 long-range scans in a tetrahedral configuration
    # m=1 short-range scan in the gap near the apex
    scans = [
        # Long-range scans (R=20)
        {'center': (0.0, 23.0, 30.0), 'radius': R_LONG, 'type': 'long'},
        {'center': (-20.0, -12.0, 30.0), 'radius': R_LONG, 'type': 'long'},
        {'center': (20.0, -12.0, 30.0), 'radius': R_LONG, 'type': 'long'},
        {'center': (0.0, 0.0, 63.0), 'radius': R_LONG, 'type': 'long'},
        # Short-range scan (R=7)
        {'center': (0.0, 0.0, 90.0), 'radius': R_SHORT, 'type': 'short'}
    ]

    # --- Verification Functions ---
    def is_contained(center, radius):
        """Checks if a sphere is fully inside the pyramid."""
        cx, cy, cz = center
        if not (cz >= radius):
            return False
        
        # Check against the four slanted faces
        max_coord_dist = HEIGHT * max(abs(cx), abs(cy))
        if max_coord_dist + (BASE_SIDE / 2) * cz > PLANE_CONSTANT - radius * PLANE_DIST_FACTOR:
            return False
        
        return True

    def check_all_scans(scan_list):
        """Verifies containment and non-overlapping for the entire list of scans."""
        for i, s1 in enumerate(scan_list):
            # Verify containment for each scan
            if not is_contained(s1['center'], s1['radius']):
                print(f"Error: Scan {i+1} at {s1['center']} is not contained.")
                return False
            
            # Verify non-overlapping with all other scans
            for j in range(i + 1, len(scan_list)):
                s2 = scan_list[j]
                dist_sq = sum([(c1 - c2)**2 for c1, c2 in zip(s1['center'], s2['center'])])
                min_dist = s1['radius'] + s2['radius']
                if dist_sq < min_dist**2:
                    print(f"Error: Scan {i+1} and Scan {j+1} are overlapping.")
                    return False
        return True

    # --- Main Logic ---
    if not check_all_scans(scans):
        print("Proposed scan configuration is invalid.")
        return

    n = sum(1 for s in scans if s['type'] == 'long')
    m = sum(1 for s in scans if s['type'] == 'short')
    
    total_scanned_volume = n * V_LONG + m * V_SHORT
    coverage_ratio = total_scanned_volume / V_PYRAMID

    # --- Output Results ---
    print("Optimal Scanning Configuration Found:")
    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print("\nCalculation of Coverage Ratio (p):")
    print(f"Volume of Pyramid = 1/3 * {BASE_SIDE}^2 * {HEIGHT} = {V_PYRAMID:.1f} m^3")
    print(f"Volume of one long-range scan = 4/3 * pi * {R_LONG}^3 = {V_LONG:.1f} m^3")
    print(f"Volume of one short-range scan = 4/3 * pi * {R_SHORT}^3 = {V_SHORT:.1f} m^3")
    print(f"Total Scanned Volume = {n} * {V_LONG:.1f} + {m} * {V_SHORT:.1f} = {total_scanned_volume:.1f} m^3")
    print(f"Coverage Ratio (p) = {total_scanned_volume:.1f} / {V_PYRAMID:.1f} = {coverage_ratio:.5f}")
    
    p_percent_str = f"{coverage_ratio:.1%}"
    print(f"\nFinal Answer (n:m:p):")
    print(f"{n}:{m}:{p_percent_str}")

solve()
<<<4:1:16.4%>>>