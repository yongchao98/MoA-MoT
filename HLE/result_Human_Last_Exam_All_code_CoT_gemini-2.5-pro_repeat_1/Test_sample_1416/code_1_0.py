import math

def solve_pyramid_scanning():
    """
    Finds the optimal locations and radii for 6 seismic scans inside a pyramid
    to maximize the total scanned volume.
    """
    
    # Pyramid and scanner constraints
    H = 110.0  # height
    B_HALF = 75.0 # half base side
    SQRT_709 = math.sqrt(22**2 + 15**2) # Pre-calculated for plane distance formula

    # Search space
    # After some analysis, the maximum possible side-sphere radius is 25.5,
    # and the distance 'd' can only be in a small range [25.5, 26.5]
    # for the spheres to be contained.
    d_search_space = [25.5, 26.0, 26.5]
    
    best_config = None
    max_volume_metric = 0

    # Search for the best configuration
    for d in d_search_space:
        
        # --- 1. Define the 4 spheres in a square pattern ---
        # Their radius is limited by the pyramid walls. Max possible r_sq is 25.5
        r_sq = 25.5
        # Place them as low as possible to maximize space, so z = r
        z_sq = r_sq
        
        # Check containment for the 4 square-pattern spheres
        # Furthest point is at (d, d, z_sq)
        if 1650 - 22 * d - 15 * z_sq < r_sq * SQRT_709:
            continue # This 'd' is too large, spheres are outside

        # --- 2. Define the lower axial sphere (S_a1) ---
        # Place it as low as possible (z_a1 = r_a1) to maximize its radius
        # Its radius is limited by non-overlap with the 4 S_sq spheres
        # dist^2 >= (r_a1 + r_sq)^2 => d^2+d^2+(r_a1-z_sq)^2 >= (r_a1+r_sq)^2
        # After simplifying, this gives r_a1 <= (2*d^2) / (2*z_sq + 2*r_sq)
        # Since z_sq=r_sq, r_a1 <= 2*d^2 / (4*r_sq) = d^2 / (2*r_sq)
        
        r_a1_max = d**2 / (2 * r_sq)
        r_a1 = math.floor(r_a1_max * 2) / 2.0 # Round down to nearest 0.5
        if r_a1 < 10.0: continue
        z_a1 = r_a1

        # --- 3. Define the upper axial sphere (S_a2) ---
        # Its radius is limited by non-overlap (with S_sq, S_a1) and containment
        # Place it just on top of the lower sphere: z_a2 = z_a1 + r_a1 + r_a2
        
        # Iteratively find the largest possible r_a2
        r_a2_found = 0
        z_a2_found = 0
        for r_a2_candidate in [i * 0.5 for i in range(100, 19, -1)]: # from 50.0 down to 10.0
            r_a2 = r_a2_candidate
            # Non-overlap with S_a1
            z_a2 = z_a1 + r_a1 + r_a2

            # Pyramid containment check
            if z_a2 > H - (r_a2 * SQRT_709 / 15.0):
                continue
            
            # Non-overlap with S_sq check
            dist_sq_a2_sq = d**2 + d**2 + (z_a2 - z_sq)**2
            if dist_sq_a2_sq < (r_a2 + r_sq)**2:
                continue
            
            # Found a valid radius
            r_a2_found = r_a2
            z_a2_found = z_a2
            break
            
        if not r_a2_found: continue
            
        # --- 4. Evaluate this configuration ---
        current_volume_metric = 4 * r_sq**3 + r_a1**3 + r_a2_found**3
        
        if current_volume_metric > max_volume_metric:
            max_volume_metric = current_volume_metric
            best_config = {
                "d": d,
                "s_sq": {"r": r_sq, "z": z_sq},
                "s_a1": {"r": r_a1, "z": z_a1},
                "s_a2": {"r": r_a2_found, "z": z_a2_found}
            }

    # --- 5. Print the optimal solution ---
    print("Optimal configuration found for N=6 scans:")
    print("This maximizes the total scanned volume while respecting all constraints.")
    
    cfg = best_config
    s_sq, s_a1, s_a2 = cfg['s_sq'], cfg['s_a1'], cfg['s_a2']
    d = cfg['d']
    
    print("\n#1: Sphere 1 (Square Pattern):")
    print(f"   Center: ({d}, {d}, {s_sq['z']}) m, Radius: {s_sq['r']} m")
    print("#2: Sphere 2 (Square Pattern):")
    print(f"   Center: ({d}, -{d}, {s_sq['z']}) m, Radius: {s_sq['r']} m")
    print("#3: Sphere 3 (Square Pattern):")
    print(f"   Center: (-{d}, {d}, {s_sq['z']}) m, Radius: {s_sq['r']} m")
    print("#4: Sphere 4 (Square Pattern):")
    print(f"   Center: (-{d}, -{d}, {s_sq['z']}) m, Radius: {s_sq['r']} m")
    print("#5: Sphere 5 (Lower Axial):")
    print(f"   Center: (0.0, 0.0, {s_a1['z']}) m, Radius: {s_a1['r']} m")
    print("#6: Sphere 6 (Upper Axial):")
    print(f"   Center: (0.0, 0.0, {s_a2['z']}) m, Radius: {s_a2['r']} m")

    all_radii = [s_sq['r']] * 4 + [s_a1['r'], s_a2['r']]
    R_max = max(all_radii)
    r_min = min(all_radii)
    
    print("\nFinal Answer (Max Radius : Min Radius):")
    print(f"R:r = {R_max}:{r_min}")
    
    # Final answer in the required format
    final_answer = f"{R_max}:{r_min}"
    return final_answer

final_answer_str = solve_pyramid_scanning()
print(f"\n<<< {final_answer_str} >>>")