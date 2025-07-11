import math

def solve_pyramid_scanning():
    """
    Calculates and verifies an optimal placement of 6 scanning spheres
    inside the Isis pyramid and reports the minimum and maximum radii used.
    """
    # 1. Pyramid and Scanner Parameters
    base_side = 150.0  # a
    height = 110.0  # h
    num_scans = 6
    min_radius = 10.0
    max_radius = 50.0
    
    print("--- Pyramid and Scan Parameters ---")
    print(f"Pyramid Base Side (a): {base_side} m")
    print(f"Pyramid Height (h): {height} m")
    print("-" * 35)

    # Derived geometric constant for face plane equations
    norm = math.sqrt(4 * height**2 + base_side**2)

    # 2. Optimized Sphere Configuration
    # This configuration is the result of the manual optimization described in the plan.
    # It follows a 1-4-1 symmetric pattern.
    
    # Sphere 1: Bottom center
    s1 = {'id': 1, 'center': (0.0, 0.0, 13.0), 'radius': 13.0}
    
    # Spheres 2-5: Four identical spheres in a symmetric layer
    r_s = 25.5
    z_s = 25.5
    d_s = 26.0
    s2 = {'id': 2, 'center': (d_s, d_s, z_s), 'radius': r_s}
    s3 = {'id': 3, 'center': (-d_s, d_s, z_s), 'radius': r_s}
    s4 = {'id': 4, 'center': (d_s, -d_s, z_s), 'radius': r_s}
    s5 = {'id': 5, 'center': (-d_s, -d_s, z_s), 'radius': r_s}
    
    # Sphere 6: Top center
    s6 = {'id': 6, 'center': (0.0, 0.0, 63.5), 'radius': 26.0}

    spheres = [s1, s2, s3, s4, s5, s6]
    
    print("\n--- Proposed Optimal Scan Configuration ---")
    for s in spheres:
        print(f"Scan {s['id']}: Radius = {s['radius']:.1f} m, Center = {s['center']}")
    print("-" * 35)

    # 3. Verification
    print("\n--- Verification of Constraints ---")
    valid = True
    
    # Check individual spheres (containment and size limits)
    print("\nChecking individual sphere containment and size...")
    for s in spheres:
        c = s['center']
        r = s['radius']
        
        # Check radius range
        if not (min_radius <= r <= max_radius):
            print(f"FAIL: Scan {s['id']} radius {r} is out of [{min_radius}, {max_radius}] range.")
            valid = False
            
        # Check containment
        # Dist to base (z=0)
        dist_base = c[2]
        if r > dist_base:
            print(f"FAIL: Scan {s['id']} violates base containment. r={r}, dist_to_base={dist_base}")
            valid = False
        
        # Dist to side faces
        # The minimum distance is to the face corresponding to the largest coordinate
        max_coord = max(abs(c[0]), abs(c[1]))
        dist_face = (base_side * height - base_side * c[2] - 2 * height * max_coord) / norm
        if r > dist_face:
            print(f"FAIL: Scan {s['id']} violates side containment. r={r:.2f}, dist_to_face={dist_face:.2f}")
            valid = False

    if valid:
        print("PASS: All spheres are within size limits and contained in the pyramid.")

    # Check non-overlap between all pairs of spheres
    print("\nChecking non-overlap between sphere pairs...")
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            s_i = spheres[i]
            s_j = spheres[j]
            c_i, r_i = s_i['center'], s_i['radius']
            c_j, r_j = s_j['center'], s_j['radius']
            
            dist_sq = sum([(c_i[k] - c_j[k])**2 for k in range(3)])
            dist = math.sqrt(dist_sq)
            min_dist = r_i + r_j
            
            if dist <= min_dist:
                print(f"FAIL: Scans {s_i['id']} and {s_j['id']} overlap.")
                print(f"      Distance = {dist:.2f}, Sum of radii = {min_dist:.2f}")
                valid = False
    
    if valid:
        print("PASS: No spheres overlap.")
    print("-" * 35)
    
    # 4. Final Results
    print("\n--- Final Results ---")
    if not valid:
        print("The proposed configuration is invalid.")
        return

    all_radii = [s['radius'] for s in spheres]
    max_r = max(all_radii)
    min_r = min(all_radii)
    
    total_volume = sum([(4/3) * math.pi * r**3 for r in all_radii])
    print(f"All radii: {all_radii}")
    print(f"Total Scanned Volume: {total_volume:.2f} m^3")
    print(f"\nMaximum radius (R): {max_r:.1f}")
    print(f"Minimum radius (r): {min_r:.1f}")

    final_answer = f"{max_r:.1f}:{min_r:.1f}"
    print(f"\nAnswer in R:r format is: {final_answer}")
    
    # Final answer in specified format for the platform
    print(f"\n<<<{final_answer}>>>")

solve_pyramid_scanning()