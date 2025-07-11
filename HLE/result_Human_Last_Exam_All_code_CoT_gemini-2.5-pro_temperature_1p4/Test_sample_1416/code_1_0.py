import numpy as np

def solve_pyramid_scanning():
    """
    Calculates the optimal placement and radii for 6 seismic scans
    inside the Isis pyramid based on a greedy 1-4-1 placement strategy.
    """
    # --- Pyramid and Scanner Parameters ---
    BASE_SIDE = 150.0  # meters
    HEIGHT = 110.0   # meters
    MIN_RADIUS = 10.0
    MAX_RADIUS = 50.0
    STEP = 0.5

    # --- Pre-calculated constants for inside-pyramid checks ---
    # Equation for slanted faces: +/- 2*h*x + a*z - a*h <= 0
    # K = sqrt((2h)^2 + a^2)
    K = np.sqrt((2 * HEIGHT)**2 + BASE_SIDE**2)

    # List to store the final sphere configurations (xc, yc, zc, r)
    spheres = []

    print("Step 1: Placing the largest possible sphere on the central axis (S1).")
    
    s1_found = False
    # Iterate downwards from max radius to find the largest possible sphere
    for r in np.arange(MAX_RADIUS, MIN_RADIUS - STEP, -STEP):
        # A sphere (0, 0, zc) with radius r is inside if:
        # 1. zc >= r
        # 2. HEIGHT * BASE_SIDE - BASE_SIDE * zc >= r * K
        # From (2), zc <= (HEIGHT * BASE_SIDE - r * K) / BASE_SIDE
        
        zc_max_by_slant = (HEIGHT * BASE_SIDE - r * K) / BASE_SIDE
        
        # Valid zc range
        zc_min = r
        zc_max = min(HEIGHT - r, zc_max_by_slant)
        
        # Find the first valid discrete position for zc
        first_valid_zc = np.ceil(zc_min / STEP) * STEP
        
        if first_valid_zc <= zc_max:
            # Found the largest sphere. Place it at the lowest possible zc to leave room.
            s1 = {'id': 'S1 (Mid)', 'xc': 0.0, 'yc': 0.0, 'zc': first_valid_zc, 'r': r}
            spheres.append(s1)
            print(f"Found S1: radius={s1['r']}m at (0.0, 0.0, {s1['zc']})m")
            s1_found = True
            break
            
    if not s1_found:
        print("Could not place the first sphere. Aborting.")
        return

    print("\nStep 2: Placing a sphere on top of S1 (S2).")
    s2_found = False
    s1 = spheres[0]
    for r in np.arange(MAX_RADIUS, MIN_RADIUS - STEP, -STEP):
        # A sphere (0, 0, zc) with radius r on top of s1 must satisfy:
        # 1. zc >= s1['zc'] + s1['r'] + r  (non-overlap)
        # 2. zc >= r (inside base)
        # 3. zc <= HEIGHT - r (inside top)
        # 4. zc <= (HEIGHT * BASE_SIDE - r * K) / BASE_SIDE (inside slant)
        
        zc_max_by_slant = (HEIGHT * BASE_SIDE - r * K) / BASE_SIDE
        
        # Valid zc range
        zc_min_non_overlap = s1['zc'] + s1['r'] + r
        zc_min_inside = r
        
        zc_min = max(zc_min_non_overlap, zc_min_inside)
        zc_max = min(HEIGHT - r, zc_max_by_slant)
        
        first_valid_zc = np.ceil(zc_min / STEP) * STEP
        
        if first_valid_zc <= zc_max:
            # Found the largest possible sphere on top.
            s2 = {'id': 'S2 (Top)', 'xc': 0.0, 'yc': 0.0, 'zc': first_valid_zc, 'r': r}
            spheres.append(s2)
            print(f"Found S2: radius={s2['r']}m at (0.0, 0.0, {s2['zc']})m")
            s2_found = True
            break

    if not s2_found:
        print("Could not place the second sphere.")

    print("\nStep 3: Placing four corner spheres (S3-S6) near the base.")
    s_corner_found = False
    s1 = spheres[0]
    for r in np.arange(MAX_RADIUS, MIN_RADIUS - STEP, -STEP):
        # We assume corner spheres touch the base for max radius, so zc = r.
        zc = r
        if zc > HEIGHT - r: continue # Basic check

        # Find if a valid placement distance 'd' exists for centers (+/-d, +/-d, zc)
        # Constraint 1: Non-overlap with S1
        # dist_sq = (d-s1['xc'])^2 + (d-s1['yc'])^2 + (zc-s1['zc'])^2 >= (r + s1['r'])^2
        # 2*d^2 + (zc-s1['zc'])^2 >= (r + s1['r'])^2
        # d^2 >= ((r + s1['r'])**2 - (zc - s1['zc'])**2) / 2
        
        # Constraint 2: Non-overlap between corner spheres (e.g. (d,d) and (d,-d))
        # dist = 2d >= 2r => d >= r
        
        # Constraint 3: Inside pyramid slanted walls
        # HEIGHT*BASE_SIDE - BASE_SIDE*zc - (2*HEIGHT)*d >= r*K
        # d <= (HEIGHT*BASE_SIDE - BASE_SIDE*zc - r*K) / (2*HEIGHT)
        
        # From C1
        c1_rhs = ((r + s1['r'])**2 - (zc - s1['zc'])**2)
        if c1_rhs < 0: continue # Should not happen in this geometry
        d_min_c1 = np.sqrt(c1_rhs / 2.0)
        
        # From C2
        d_min_c2 = r

        d_min = max(d_min_c1, d_min_c2)
        
        # From C3
        d_max_c3_num = HEIGHT * BASE_SIDE - BASE_SIDE * zc - r * K
        if d_max_c3_num < 0: continue
        d_max = d_max_c3_num / (2 * HEIGHT)
        
        first_valid_d = np.ceil(d_min / STEP) * STEP
        
        if first_valid_d <= d_max:
            d = first_valid_d
            print(f"Found corner spheres: radius={r}m")
            print(f"  - Center z-coordinate: {zc}m")
            print(f"  - Center xy-distance 'd': {d}m")
            
            centers = [(d, d), (d, -d), (-d, d), (-d, -d)]
            ids = ['S3', 'S4', 'S5', 'S6']
            for i in range(4):
                corner_sphere = {'id': f'{ids[i]} (Corner)', 'xc': centers[i][0], 'yc': centers[i][1], 'zc': zc, 'r': r}
                spheres.append(corner_sphere)
            s_corner_found = True
            break

    if not s_corner_found:
        print("Could not place the four corner spheres with this strategy.")

    print("\n--- Final Scanning Configuration (N=6) ---")
    if len(spheres) == 6:
        all_radii = []
        for s in spheres:
            print(f"Scan {s['id']:<12}: Center=({s['xc']:.1f}, {s['yc']:.1f}, {s['zc']:.1f}), Radius={s['r']:.1f}")
            all_radii.append(s['r'])
        
        max_r = max(all_radii)
        min_r = min(all_radii)
        
        print("\n--- Result ---")
        print(f"The calculated radii are: {all_radii}")
        print(f"Maximum radius (R) = {max_r:.1f} m")
        print(f"Minimum radius (r) = {min_r:.1f} m")
        print(f"Final Answer (R:r) = {max_r:.1f}:{min_r:.1f}")
        
        # Required final answer format
        print(f"\n<<<{max_r:.1f}:{min_r:.1f}>>>")

    else:
        print(f"\nCould not find a valid configuration for N=6. Found {len(spheres)} spheres.")

# Run the calculation
solve_pyramid_scanning()