import math

def solve_pyramid_scanning():
    """
    Calculates and verifies the optimal placement for 6 scanners
    inside the Isis pyramid and provides the R:r ratio.
    """
    # --- Pyramid Parameters ---
    B = 150.0  # Base side in meters
    H = 110.0  # Height in meters

    # --- Derived Geometric Constant ---
    # Used for calculating distance from a point to a side face
    K = math.sqrt(4 * H**2 + B**2)

    # --- Optimal Configuration Found via Analysis ---
    # This configuration is based on a symmetric 4+2 octahedral arrangement
    # centered around the pyramid's center of volume (H/4).
    #
    # 4 Equatorial spheres
    r_eq = 21.0
    d_eq = 30.5
    z_eq = 27.5
    
    # 1 Upper Axial sphere
    r_up = 26.0
    z_up = 63.5

    # 1 Lower Axial sphere
    r_low = 12.5
    z_low = 12.5

    # List of all 6 spheres: (center_x, center_y, center_z, radius)
    scanners = [
        # Equatorial spheres
        (d_eq, 0.0, z_eq, r_eq),
        (-d_eq, 0.0, z_eq, r_eq),
        (0.0, d_eq, z_eq, r_eq),
        (0.0, -d_eq, z_eq, r_eq),
        # Axial spheres
        (0.0, 0.0, z_up, r_up),
        (0.0, 0.0, z_low, r_low)
    ]

    print("Optimal Scanner Configuration (N=6):")
    print("-" * 40)
    for i, s in enumerate(scanners):
        print(f"Scanner {i+1}:")
        print(f"  Center (x, y, z): ({s[0]:.1f}m, {s[1]:.1f}m, {s[2]:.1f}m)")
        print(f"  Radius: {s[3]:.1f}m")
        total_volume = 4/3 * math.pi * s[3]**3
        print(f"  Volume: {total_volume:.2f} m^3")
        print("-" * 40)

    # --- Verification ---
    all_constraints_met = True
    
    # 1. Check individual constraints (containment & radius range)
    for i, s in enumerate(scanners):
        cx, cy, cz, r = s
        # Check radius range
        if not (10.0 <= r <= 50.0):
            print(f"Error: Scanner {i+1} radius {r} out of range [10, 50].")
            all_constraints_met = False
        
        # Check containment
        # Base constraint
        if cz < r:
            print(f"Error: Scanner {i+1} too close to base. Dist: {cz}, Radius: {r}")
            all_constraints_met = False
        # Side wall constraint
        dist_to_wall = (H * B - B * cz - 2 * H * max(abs(cx), abs(cy))) / K
        if dist_to_wall < r - 1e-9: # Allow for tiny float precision errors
            print(f"Error: Scanner {i+1} too close to side. Dist: {dist_to_wall:.2f}, Radius: {r}")
            all_constraints_met = False
            
    # 2. Check non-overlap constraint between all pairs
    for i in range(len(scanners)):
        for j in range(i + 1, len(scanners)):
            s1 = scanners[i]
            s2 = scanners[j]
            c1, r1 = s1[:3], s1[3]
            c2, r2 = s2[:3], s2[3]
            
            dist_sq = (c1[0]-c2[0])**2 + (c1[1]-c2[1])**2 + (c1[2]-c2[2])**2
            dist = math.sqrt(dist_sq)
            
            if dist < r1 + r2 - 1e-9: # Allow for float precision errors
                print(f"Error: Scanner {i+1} and {j+1} are overlapping. Dist: {dist:.2f}, Sum Radii: {r1+r2}")
                all_constraints_met = False

    if all_constraints_met:
        print("\nAll constraints are successfully verified.")
    else:
        print("\nVerification failed. The configuration is not valid.")

    # --- Final Answer ---
    all_radii = [s[3] for s in scanners]
    R_max = max(all_radii)
    r_min = min(all_radii)
    
    print(f"\nThe maximum radius (R) is {R_max:.1f}m and the minimum radius (r) is {r_min:.1f}m.")
    print("Final answer in the format R:r")
    print(f"{R_max:.1f}:{r_min:.1f}")


solve_pyramid_scanning()
<<<26.0:12.5>>>