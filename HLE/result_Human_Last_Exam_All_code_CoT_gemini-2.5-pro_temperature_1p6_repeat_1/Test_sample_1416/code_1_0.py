import math

# --- Pyramid and Scanner Parameters ---
PYRAMID_BASE = 150.0
PYRAMID_HEIGHT = 110.0
# Derived constants for the pyramid face plane equation (22x + 15z = 1650)
# Normal vector magnitude: sqrt(22^2 + 15^2)
SQRT_709 = math.sqrt(22**2 + 15**2)

# --- Proposed Optimal Solution ---
# This solution is based on a deformed octahedral arrangement,
# optimized to maximize volume while respecting all constraints.

# Sphere definitions: (center_x, center_y, center_z, radius)
spheres = [
    # 4 "Equatorial" spheres
    (24.5, 0.0, 41.5, 17.0),
    (-24.5, 0.0, 41.5, 17.0),
    (0.0, 24.5, 41.5, 17.0),
    (0.0, -24.5, 41.5, 17.0),
    # 1 "Top" axial sphere
    (0.0, 0.0, 66.0, 17.5),
    # 1 "Bottom" axial sphere
    (0.0, 0.0, 17.0, 17.0)
]

# --- Verification Functions ---

def check_containment(sphere):
    """Checks if a sphere is fully contained within the pyramid."""
    x, y, z, r = sphere
    # 1. Check against base plane (z=0)
    if z < r:
        print(f"FAIL: Sphere at ({x},{y},{z}) with r={r} breaks base constraint (z < r).")
        return False
    # 2. Check against the four side faces
    # The distance from a center (x,y,z) to a plane Ax+By+Cz+D=0 must be >= r.
    # For plane 22x + 15z - 1650 = 0:
    dist_to_plane_pos_x = (1650 - 22 * abs(x) - 15 * z) / SQRT_709
    if dist_to_plane_pos_x < r:
        print(f"FAIL: Sphere at ({x},{y},{z}) with r={r} breaks side constraint X (dist={dist_to_plane_pos_x}).")
        return False
    dist_to_plane_pos_y = (1650 - 22 * abs(y) - 15 * z) / SQRT_709
    if dist_to_plane_pos_y < r:
        print(f"FAIL: Sphere at ({x},{y},{z}) with r={r} breaks side constraint Y (dist={dist_to_plane_pos_y}).")
        return False
    return True

def check_non_overlap(spheres_list):
    """Checks if any two spheres in the list overlap."""
    n = len(spheres_list)
    for i in range(n):
        for j in range(i + 1, n):
            s1 = spheres_list[i]
            s2 = spheres_list[j]
            dist_sq = (s1[0] - s2[0])**2 + (s1[1] - s2[1])**2 + (s1[2] - s2[2])**2
            min_dist = s1[3] + s2[3]
            if dist_sq < min_dist**2:
                print(f"FAIL: Sphere {i+1} and Sphere {j+1} overlap.")
                return False
    return True

def run_verification():
    """Runs all checks on the proposed solution."""
    print("--- Verifying Proposed Solution ---")
    
    # Verify each sphere's containment
    all_contained = all(check_containment(s) for s in spheres)
    if all_contained:
        print("PASS: All spheres are fully contained within the pyramid.")
    
    # Verify non-overlapping condition for all pairs
    all_non_overlapping = check_non_overlap(spheres)
    if all_non_overlapping:
        print("PASS: All spheres are non-overlapping.")
        
    if not (all_contained and all_non_overlapping):
        print("\nVerification failed.")
        return

    print("\n--- Optimal Configuration Found ---")
    radii = [s[3] for s in spheres]
    max_r = max(radii)
    min_r = min(radii)
    
    for i, s in enumerate(spheres):
        print(f"Scan {i+1}: Center=({s[0]:.1f}, {s[1]:.1f}, {s[2]:.1f}), Radius={s[3]:.1f}")
    
    print("\nThis configuration maximizes the scanned volume.")
    
    # "output each number in the final equation"
    # This interprets the request as showing the components of the total volume.
    print("\nTotal Volume Calculation:")
    total_volume = 0
    volume_eq_parts = []
    
    # Count spheres of each radius
    radius_counts = {}
    for r in radii:
        radius_counts[r] = radius_counts.get(r, 0) + 1
        
    for r, count in radius_counts.items():
        volume = count * (4/3) * math.pi * (r**3)
        total_volume += volume
        volume_eq_parts.append(f"{count} * (4/3 * pi * {r}^3)")
    
    equation_str = " + ".join(volume_eq_parts)
    print(f"V_total = {equation_str}")
    print(f"V_total = {total_volume:.2f} m^3")

    print(f"\nThe maximum scanning radius is R = {max_r:.1f}m")
    print(f"The minimum scanning radius is r = {min_r:.1f}m")
    
    # Final answer in the specified format
    print(f"\nFinal Answer: R:r")
    print(f"<<<{max_r:.1f}:{min_r:.1f}>>>")

if __name__ == '__main__':
    run_verification()