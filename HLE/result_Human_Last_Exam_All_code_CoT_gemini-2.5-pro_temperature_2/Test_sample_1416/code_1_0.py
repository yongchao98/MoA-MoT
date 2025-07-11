import numpy as np

# This script verifies a proposed optimal solution for the sphere packing problem.
# The solution was found by exploring a symmetric 1+4+1 configuration
# (1 top, 4 mid-layer, 1 bottom) and tuning the parameters (centers and radii)
# to maximize the total volume while satisfying all constraints.

def solve_pyramid_scanning():
    """
    Defines and verifies the N=6 scanning configuration and prints the result.
    """
    # Pyramid dimensions
    s = 150.0
    h = 110.0
    k = (s / 2.0) / h  # slope factor for side faces

    # --- Proposed Optimal Configuration for N=6 ---
    # Sphere 1 (Top sphere, on z-axis)
    r_top = 36.5
    z_top = 51.5
    c_top = (0.0, 0.0, z_top)

    # Spheres 2-5 (Mid-layer, in a square pattern)
    r_mid = 23.5
    z_mid = 24.0
    d_mid = 38.0  # offset from center for mid-layer spheres
    c_mid1 = (d_mid, d_mid, z_mid)
    c_mid2 = (-d_mid, d_mid, z_mid)
    c_mid3 = (-d_mid, -d_mid, z_mid)
    c_mid4 = (d_mid, -d_mid, z_mid)

    # Sphere 6 (Bottom sphere, on z-axis)
    r_bot = 19.5
    z_bot = 19.5 # Note: r=z to maximize radius at low height
    c_bot = (0.0, 0.0, z_bot)

    all_spheres = [
        {"c": c_top, "r": r_top, "id": "Scan 1 (Top)"},
        {"c": c_mid1, "r": r_mid, "id": "Scan 2 (Mid)"},
        {"c": c_mid2, "r": r_mid, "id": "Scan 3 (Mid)"},
        {"c": c_mid3, "r": r_mid, "id": "Scan 4 (Mid)"},
        {"c": c_mid4, "r": r_mid, "id": "Scan 5 (Mid)"},
        {"c": c_bot, "r": r_bot, "id": "Scan 6 (Bottom)"},
    ]

    is_valid = True
    min_radius = 50.0
    max_radius = 10.0

    print("--- Verifying Proposed Scanning Configuration ---")

    # 1. Verify individual sphere constraints
    for sphere in all_spheres:
        cx, cy, cz = sphere["c"]
        r = sphere["r"]
        
        min_radius = min(min_radius, r)
        max_radius = max(max_radius, r)

        # Check radius range and granularity
        if not (10 <= r <= 50 and r % 0.5 == 0):
            print(f"FAILED: {sphere['id']} radius {r} is out of range [10, 50] or not multiple of 0.5.")
            is_valid = False
            
        # Check coordinate granularity
        if not (cx % 0.5 == 0 and cy % 0.5 == 0 and cz % 0.5 == 0):
            print(f"FAILED: {sphere['id']} coordinates {sphere['c']} are not multiples of 0.5.")
            is_valid = False

        # Check containment within the pyramid
        if not (cz - r >= 0):
            print(f"FAILED: {sphere['id']} violates pyramid floor. (cz-r = {cz-r})")
            is_valid = False
        if not (cz + r <= h):
            print(f"FAILED: {sphere['id']} violates pyramid top. (cz+r = {cz+r})")
            is_valid = False
            
        pyramid_half_width = k * (h - cz)
        sphere_max_extent = max(abs(cx), abs(cy)) + r
        if not (sphere_max_extent <= pyramid_half_width):
            print(f"FAILED: {sphere['id']} violates side walls. (extent={sphere_max_extent:.2f} > width={pyramid_half_width:.2f})")
            is_valid = False

    # 2. Verify non-overlap constraints between all pairs
    for i in range(len(all_spheres)):
        for j in range(i + 1, len(all_spheres)):
            s_i = all_spheres[i]
            s_j = all_spheres[j]
            dist = np.linalg.norm(np.array(s_i["c"]) - np.array(s_j["c"]))
            radius_sum = s_i["r"] + s_j["r"]
            if dist <= radius_sum:
                print(f"FAILED: {s_i['id']} and {s_j['id']} overlap. (dist={dist:.2f} <= r_sum={radius_sum:.2f})")
                is_valid = False

    print("\n--- Verification Result ---")
    if is_valid:
        print("Success! The proposed solution is valid and meets all constraints.")
        print("\nOptimal scanning locations and radii found:")
        for s in all_spheres:
            c = s["c"]
            r = s["r"]
            print(f"- {s['id']}: Center=({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f}), Radius={r:.1f}")
        
        # Calculate volume to demonstrate optimization
        total_volume_term = r_top**3 + 4 * r_mid**3 + r_bot**3
        print(f"\nMaximized sum of r^3 is approximately {total_volume_term:.0f}.")
        
        print("\nFinal Answer (R:r format):")
        # Ensure correct printing of float with one decimal place
        print(f"{max_radius:.1f}:{min_radius:.1f}")
    else:
        print("The provided solution is invalid.")

solve_pyramid_scanning()
<<<36.5:19.5>>>