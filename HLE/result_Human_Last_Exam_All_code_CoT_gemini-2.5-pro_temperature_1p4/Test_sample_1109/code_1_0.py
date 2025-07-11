import math

def solve_seismic_scanning():
    """
    Calculates the optimal number of long and short range scans and the
    resulting coverage of the Isis pyramid.
    """

    # 1. Define constants and geometric properties
    S_BASE = 150.0
    HEIGHT = 110.0
    R_LONG = 20.0
    R_SHORT = 7.0
    STEP = 0.5
    
    # Pre-calculated constant for the pyramid side plane equation
    # Plane: 22x + 15z - 1650 = 0. Normal vector component ratio is 22:15.
    SQRT_709 = math.sqrt(22**2 + 15**2)

    # 2. Define the optimal long-range scan configuration based on analysis
    long_range_placements = [
        {'center': (20.0, 20.0, 20.0), 'radius': R_LONG},
        {'center': (-20.0, 20.0, 20.0), 'radius': R_LONG},
        {'center': (20.0, -20.0, 20.0), 'radius': R_LONG},
        {'center': (-20.0, -20.0, 20.0), 'radius': R_LONG},
        {'center': (0.0, 0.0, 60.0), 'radius': R_LONG}
    ]
    
    # 3. Pack the short-range spheres greedily
    short_range_placements = []
    r = R_SHORT
    all_placed_spheres = long_range_placements + short_range_placements
    
    # Iterate through possible z-coordinates from bottom up
    z_coords = [i * STEP for i in range(int(r / STEP), int(HEIGHT / STEP))]

    for cz in z_coords:
        # Determine max |cx|, |cy| for a center at this height to be inside
        max_c_abs = (1650 - 15 * cz - r * SQRT_709) / 22
        if max_c_abs < 0:
            continue
        
        max_dist_steps = int(max_c_abs / STEP)
        
        # Generate candidate (cx, cy) points for this z-level, ordered center-out
        points_to_check = set()
        for dist_step in range(max_dist_steps + 1):
            dist = dist_step * STEP
            for c_coord_step in range(-dist_step, dist_step + 1):
                c_coord = c_coord_step * STEP
                if dist_step == abs(c_coord_step): # only check corners once
                    points_to_check.add((c_coord, dist))
                    points_to_check.add((c_coord, -dist))
                else:
                    points_to_check.add((c_coord, dist))
                    points_to_check.add((c_coord, -dist))
                    points_to_check.add((dist, c_coord))
                    points_to_check.add((-dist, c_coord))

        sorted_points = sorted(list(points_to_check), key=lambda p: p[0]**2 + p[1]**2)

        # Filter spheres that are too far on z-axis to possibly overlap
        relevant_spheres = [s for s in all_placed_spheres if abs(s['center'][2] - cz) < s['radius'] + r]

        for cx, cy in sorted_points:
            candidate_center = (cx, cy, cz)
            is_overlapped = False
            for s in relevant_spheres:
                dist_sq = ( (candidate_center[0] - s['center'][0])**2 + 
                            (candidate_center[1] - s['center'][1])**2 + 
                            (candidate_center[2] - s['center'][2])**2 )
                if dist_sq < (r + s['radius'])**2:
                    is_overlapped = True
                    break
            
            if not is_overlapped:
                new_sphere = {'center': candidate_center, 'radius': r}
                short_range_placements.append(new_sphere)
                all_placed_spheres.append(new_sphere)
                relevant_spheres.append(new_sphere)

    # 4. Calculate final results
    n = len(long_range_placements)
    m = len(short_range_placements)

    v_pyramid = (S_BASE**2 * HEIGHT) / 3.0
    v_long_sphere = (4.0 / 3.0) * math.pi * (R_LONG**3)
    v_short_sphere = (4.0 / 3.0) * math.pi * (R_SHORT**3)
    v_scanned = n * v_long_sphere + m * v_short_sphere

    coverage_ratio = (v_scanned / v_pyramid) * 100.0

    # 5. Output the numbers and the final formatted answer
    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print(f"Total coverage ratio (p): {coverage_ratio:.1f}%")
    
    final_answer = f"{n}:{m}:{coverage_ratio:.1f}"
    print("\nFinal Answer Format:")
    print(final_answer)
    return final_answer

final_result_str = solve_seismic_scanning()
print(f"<<<{final_result_str}>>>")