import math

def solve_pyramid_scanning():
    """
    Calculates the optimal placement and radii for 6 seismic scanners
    in a square pyramid using a 4+1+1 symmetric arrangement.
    """
    # Pyramid dimensions
    base_side = 150.0  # b
    height = 110.0     # h

    # Pre-calculate a constant used for distance to side planes
    # The equation for a side plane is 2*|x|/b + z/h = 1
    # The distance to this plane from a point (cx,cy,cz) is
    # (1 - 2*max(|cx|,|cy|)/b - cz/h) / D_val
    D_val = math.sqrt(4 / base_side**2 + 1 / height**2)

    def get_max_radius(x, y, z):
        """Calculates the max radius of a sphere centered at (x,y,z)
           that fits inside the pyramid."""
        if z < 0:
            return 0
        # The point is inside if 2*max(|x|,|y|)/b + z/h <= 1
        if 2 * max(abs(x), abs(y)) / base_side + z / height > 1:
            return 0
        
        # Distance to the base plane (z=0)
        dist_to_base = z
        # Distance to the closest side plane
        dist_to_side = (1 - 2 * max(abs(x), abs(y)) / base_side - z / height) / D_val
        
        return min(dist_to_base, dist_to_side)

    # --- Step 1: Optimize the 4 spheres in a square plane ---
    # We place 4 spheres at (+/-d, +/-d, d) with radius r=d.
    # This configuration is tangent to the base and side walls simultaneously.
    # The limiting constraint is d*D_val <= 1 - 2*d/b - d/h
    # d * (D_val + 2/b + 1/h) <= 1
    d_max = 1 / (D_val + 2/base_side + 1/height)
    # Round down to the nearest 0.5m
    d = math.floor(d_max * 2) / 2
    r4 = d
    
    four_spheres_centers = [
        (d, d, d), (-d, d, d), (d, -d, d), (-d, -d, d)
    ]
    four_spheres_radii = [r4] * 4

    # --- Step 2: Optimize the bottom sphere on the central axis ---
    # Place a sphere at c=(0,0,z_b) with radius r_b.
    # It must not overlap with the 4 spheres. Its radius is also limited by the pyramid floor (r_b <= z_b).
    # To maximize volume, we assume r_b = z_b.
    # Non-overlap: dist((0,0,z_b), (d,d,d)) >= r_b + r4
    # We solve for z_b = r_b where it's tangent to the 4 spheres:
    # (z_b + r4)^2 = dist^2 = (d-0)^2 + (d-0)^2 + (d-z_b)^2
    # z_b^2 + 2*z_b*d + d^2 = 2*d^2 + d^2 - 2*z_b*d + z_b^2
    # 4 * z_b * d = 2 * d^2  => z_b = d / 2
    z_b = d / 2
    # Round radius down to nearest 0.5
    r_b = math.floor(z_b * 2) / 2
    # Adjust center z to match radius, as r_b=z_b is the tightest constraint
    c_b = (0, 0, r_b)

    # --- Step 3: Optimize the top sphere on the central axis ---
    # Place a sphere at c=(0,0,z_t) with radius r_t.
    # It must not overlap with the 4 spheres or the bottom sphere.
    # We search for the optimal z_t that maximizes r_t.
    best_rt = 0
    best_zt = 0
    # Search for z_t in a reasonable range above the 4 spheres
    for z_t_candidate in [i * 0.5 for i in range(int(d * 2), int(height * 2))]:
        # Containment constraint from pyramid walls
        r_contain = get_max_radius(0, 0, z_t_candidate)
        
        # Non-overlap constraint with the 4 spheres
        dist_to_4_spheres = math.sqrt(2 * d**2 + (z_t_candidate - d)**2)
        r_overlap_4 = dist_to_4_spheres - r4
        
        # Non-overlap constraint with bottom sphere
        dist_to_bottom = z_t_candidate - c_b[2]
        r_overlap_b = dist_to_bottom - r_b

        current_rt = min(r_contain, r_overlap_4, r_overlap_b)
        
        if current_rt > best_rt:
            best_rt = current_rt
            best_zt = z_t_candidate
            
    # Round radius down to nearest 0.5
    r_t = math.floor(best_rt * 2) / 2
    c_t = (0, 0, best_zt)

    # --- Step 4: Consolidate and print results ---
    all_centers = four_spheres_centers + [c_t, c_b]
    all_radii = four_spheres_radii + [r_t, r_b]
    
    print("Optimal scanning locations and radii for N=6:")
    print("-" * 50)
    for i, (center, radius) in enumerate(zip(all_centers, all_radii)):
        print(f"Scan {i+1}: Center=({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}) m, Radius={radius:.1f} m")
    print("-" * 50)

    max_radius = max(all_radii)
    min_radius = min(all_radii)
    
    print(f"Maximum radius (R): {max_radius:.1f} m")
    print(f"Minimum radius (r): {min_radius:.1f} m")
    print("\nFinal Answer Format (R:r)")
    print(f"{max_radius:.1f}:{min_radius:.1f}")
    
    # Final answer for the submission system
    final_answer_str = f"<<<{max_radius:.1f}:{min_radius:.1f}>>>"
    return final_answer_str

# Execute the function and capture the final answer string
final_answer = solve_pyramid_scanning()
print(final_answer)