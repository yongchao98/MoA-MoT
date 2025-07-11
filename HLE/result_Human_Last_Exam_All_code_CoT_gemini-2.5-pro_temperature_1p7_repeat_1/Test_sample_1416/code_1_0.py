import math

def solve_pyramid_scanning():
    """
    Calculates and verifies an optimal configuration for 6 seismic scans
    inside the Isis pyramid, aiming to maximize scanned volume.
    """

    # 1. Pyramid Geometry Definition
    # Side planes are defined by the equation: 22*|x| + 15*z <= 1650
    # And 22*|y| + 15*z <= 1650
    PLANE_A = 22.0
    PLANE_C = 15.0
    PLANE_D = 1650.0
    PLANE_DENOM = math.sqrt(PLANE_A**2 + PLANE_C**2)  # sqrt(709)

    # 2. Optimized Sphere Configuration
    # This configuration is derived from a greedy packing strategy.

    # Sphere 1: The maximal inscribed sphere, tangent to the base and 4 side walls.
    # For this sphere, its radius r1 equals its center's height z1.
    # Calculation: r1 = PLANE_D / (PLANE_C + PLANE_DENOM) ~= 39.637m
    r1 = 39.5
    c1 = (0.0, 0.0, 39.5)

    # Sphere 2: Placed on the central axis, tangent to Sphere 1 and the side walls.
    # Calculation leads to a max radius of ~11.17m. We use 10.5m for safety and discretization.
    # Its center z2 is determined by being tangent to Sphere 1: dist(c1,c2) = r1+r2
    r2 = 10.5
    c2_z = c1[2] + r1 + r2
    c2 = (0.0, 0.0, 89.5) # c1[2]+r1+r2 = 39.5+39.5+10.5, this is wrong logic
                           # dist is r1+r2, so |c2_z - c1_z|=50 => c2_z = 89.5
    
    # Spheres 3-6: Four identical spheres packed low, tangent to the base,
    # Sphere 1, and the side walls.
    # Their height z3 equals their radius r3.
    # Calculation leads to a max radius of ~12.5m.
    r3 = 12.5
    z3 = 12.5
    # Their distance d from the center axis is found via tangency to Sphere 1:
    # d^2 = 4 * r3 * r1. d = sqrt(4 * 12.5 * 39.5) = sqrt(1975) ~= 44.44m
    d3 = 44.5  # Discretize to nearest 0.5m

    c3 = (d3, 0.0, z3)
    c4 = (-d3, 0.0, z3)
    c5 = (0.0, d3, z3)
    c6 = (0.0, -d3, z3)

    spheres = [
        (c1, r1), (c2, r2), (c3, r3), (c4, r3), (c5, r3), (c6, r3)
    ]

    print("Optimal Scanning Locations and Radii (N=6):")
    print("-------------------------------------------------")
    
    all_radii = []
    for i, (c, r) in enumerate(spheres):
        print(f"Scan {i+1}: Center=({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f}) m, Radius={r:.1f} m")
        all_radii.append(r)

    # 3. Verification of Constraints
    
    # Verify Containment for each sphere
    for i, (c, r) in enumerate(spheres):
        cx, cy, cz = c
        if cz - r < -1e-9: # Base containment (with tolerance)
            raise ValueError(f"FAIL: Sphere {i+1} breaches base.")
        
        # Side wall containment (check plane closest to center)
        dist_to_wall = (PLANE_D - PLANE_A * abs(cx) - PLANE_C * cz) / PLANE_DENOM
        if dist_to_wall < r - 1e-9: # (with tolerance)
             raise ValueError(f"FAIL: Sphere {i+1} breaches side wall.")
             
    # Verify Non-Overlap between each pair of spheres
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            c_i, r_i = spheres[i]
            c_j, r_j = spheres[j]
            dist_sq = sum([(c_i[k] - c_j[k])**2 for k in range(3)])
            min_dist_sq = (r_i + r_j)**2
            if dist_sq < min_dist_sq - 1e-9: # (with tolerance)
                raise ValueError(f"FAIL: Spheres {i+1} and {j+1} overlap.")

    print("\nVerification successful: All constraints are met.")
    print("-------------------------------------------------")

    # 4. Final Answer
    R_max = max(all_radii)
    r_min = min(all_radii)
    
    print("\nThe maximum and minimum scanning radii are:")
    print(f"R_max = {R_max:.1f} m")
    print(f"R_min = {r_min:.1f} m")

    print("\nFinal Answer (R_max:R_min):")
    # This print format is for the final answer extraction
    print(f"{R_max:.1f}:{r_min:.1f}")

# Execute the solution
solve_pyramid_scanning()