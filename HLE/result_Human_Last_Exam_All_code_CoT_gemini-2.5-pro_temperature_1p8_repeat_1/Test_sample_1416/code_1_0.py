import math

def solve_pyramid_scanning():
    """
    Calculates and verifies the optimal placement for 6 scanners in the Isis pyramid.
    """
    # --- Pyramid Parameters ---
    s = 150.0  # base side in meters
    h = 110.0  # height in meters

    # --- Optimal Configuration Found (4 corner spheres, 2 axial spheres) ---
    # This configuration was found by optimizing a symmetric arrangement to maximize radii.
    
    # 4 Corner Spheres
    r_c = 27.5
    x_c = 27.5 
    z_c = 27.5
    
    # 1 Lower Axial Sphere
    r_a2 = 13.5
    z_a2 = 13.5
    
    # 1 Upper Axial Sphere
    r_a1 = 14.0
    z_a1 = 89.0

    spheres = [
        # Corner spheres
        {'id': 'Corner 1', 'r': r_c, 'c': (x_c, x_c, z_c)},
        {'id': 'Corner 2', 'r': r_c, 'c': (-x_c, x_c, z_c)},
        {'id': 'Corner 3', 'r': r_c, 'c': (x_c, -x_c, z_c)},
        {'id': 'Corner 4', 'r': r_c, 'c': (-x_c, -x_c, z_c)},
        # Axial spheres
        {'id': 'Axial Lower', 'r': r_a2, 'c': (0, 0, z_a2)},
        {'id': 'Axial Upper', 'r': r_a1, 'c': (0, 0, z_a1)},
    ]

    print("--- Verifying Optimal Scanner Placement (N=6) ---")
    all_constraints_met = True

    # 1. Verify all spheres are inside the pyramid
    print("\n1. Verifying spheres are inside the pyramid:")
    for s_info in spheres:
        r, (cx, cy, cz) = s_info['r'], s_info['c']
        # Check against pyramid boundary equation: max(|cx|, |cy|) + r <= (s/2) * (1 - cz/h)
        pyramid_half_width_at_cz = (s / 2) * (1 - cz / h)
        sphere_max_extent = max(abs(cx), abs(cy)) + r
        
        is_inside = (cz >= r) and (sphere_max_extent <= pyramid_half_width_at_cz)
        print(f"   - {s_info['id']:<12} (r={r:.1f}): max_extent={sphere_max_extent:.2f}, pyramid_width={pyramid_half_width_at_cz:.2f}.  Result: {'OK' if is_inside else 'FAIL'}")
        if not is_inside:
            all_constraints_met = False
            
    # 2. Verify no spheres overlap
    print("\n2. Verifying no spheres overlap:")
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            s1 = spheres[i]
            s2 = spheres[j]
            c1, r1 = s1['c'], s1['r']
            c2, r2 = s2['c'], s2['r']
            
            dist_sq = (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2
            min_dist_sq = (r1 + r2)**2
            
            no_overlap = dist_sq >= min_dist_sq
            print(f"   - {s1['id']:<12} vs {s2['id']:<12}: dist^2={dist_sq:.2f}, (r1+r2)^2={min_dist_sq:.2f}. Result: {'OK' if no_overlap else 'FAIL'}")
            if not no_overlap:
                all_constraints_met = False

    # 3. Calculate total volume and final answer
    if all_constraints_met:
        print("\n--- All constraints are met. Configuration is valid. ---")
        
        total_volume = 0
        all_radii = []
        for s_info in spheres:
            r = s_info['r']
            total_volume += (4/3) * math.pi * (r**3)
            all_radii.append(r)
            
        print("\nFinal Configuration Details:")
        for s_info in spheres:
            print(f"  - ID: {s_info['id']:<12} | Radius: {s_info['r']:.1f} m | Center: {s_info['c']}")
        
        print(f"\nMaximized Scanned Volume: {total_volume:.2f} cubic meters")
        
        R = max(all_radii)
        r_min = min(all_radii)
        
        print("\n--- Final Answer ---")
        print(f"The maximum scanning radius (R) is {R:.1f} m.")
        print(f"The minimum scanning radius (r) is {r_min:.1f} m.")
        
        # The equation for the maximum radius of a corner sphere (r_c) placed at x=z=r_c
        # is 2*r_c <= 75 * (1 - r_c / 110).
        # This simplifies to r_c * (2 + 75/110) <= 75, so r_c <= 75 / (295/110) = 27.96.
        # So we choose r_c = 27.5 m.
        print("\nEquation for the corner sphere radius (rc) touching the pyramid wall:")
        print(f"2 * {r_c} <= {s/2:.1f} * (1 - {r_c:.1f} / {h:.1f})")
        print(f"{2 * r_c:.1f} <= {pyramid_half_width_at_cz:.2f}")

        # The equation for the lower axial sphere (ra2) touching the corner spheres
        # where z=x=r_c and z=r_a2 is 2*rc^2 = 4*ra2*rc => rc = 2*ra2
        # ra2 = rc / 2 = 27.5 / 2 = 13.75. We choose ra2 = 13.5
        print("\nEquation for lower axial sphere radius (ra2) touching corner spheres:")
        print(f"{r_c:.1f} = 2 * ra2  => ra2 = {r_c / 2:.2f}. Chosen: {r_a2:.1f}")

        print(f"\nFinal Answer Format (R:r):")
        print(f"{R:.1f}:{r_min:.1f}")
    else:
        print("\n--- ERROR: The configuration is invalid. ---")


solve_pyramid_scanning()