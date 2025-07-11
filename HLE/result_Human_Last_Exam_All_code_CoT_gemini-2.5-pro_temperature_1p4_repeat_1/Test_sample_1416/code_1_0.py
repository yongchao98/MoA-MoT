import numpy as np

def solve_pyramid_scanning():
    """
    Finds optimal locations and radii for 6 scans in a pyramid.
    """
    # 1. Define Geometry and Constraints
    H = 110.0  # Height of the pyramid
    A = 150.0  # Base side length
    B_HALF = A / 2.0  # Half base side
    MIN_R, MAX_R = 10.0, 50.0
    STEP = 0.5

    # Norm of the normal vector to the side faces
    NORM_N = np.sqrt(H**2 + B_HALF**2)

    def floor_to_step(val, step):
        """Rounds a value down to the nearest step."""
        return np.floor(val / step) * step

    def max_radius_in_pyramid(c):
        """Calculates max radius for a sphere at center c to be in the pyramid."""
        cx, cy, cz = c
        # The pyramid is defined by 5 planes:
        # z >= 0
        # H*x + B_HALF*z - H*B_HALF <= 0, etc.
        d_base = cz
        d_faces = (H * B_HALF - H * abs(cx) - B_HALF * cz) / NORM_N
        
        radius = min(d_base, d_faces)
        return radius if radius > 0 else 0

    spheres = []

    # 2. Find Sphere 1 (Largest inscribed sphere on the central axis)
    best_r1, best_c1 = 0, None
    for z_coord in np.arange(STEP, H, STEP):
        c = np.array([0, 0, z_coord])
        r = max_radius_in_pyramid(c)
        if r > best_r1:
            best_r1, best_c1 = r, c
    
    r1 = floor_to_step(best_r1, STEP)
    c1 = np.round(best_c1 / STEP) * STEP
    spheres.append({'c': c1, 'r': r1})

    # 3. Find Spheres 2-5 (Symmetric corner spheres)
    best_r2, best_c2_quadrant = 0, None
    # Search for the best sphere in the +x, +y quadrant
    for x_coord in np.arange(STEP, B_HALF, STEP):
        for z_coord in np.arange(STEP, H, STEP):
            c_cand = np.array([x_coord, x_coord, z_coord])
            
            # Constraint: Pyramid wall
            r_geom = max_radius_in_pyramid(c_cand)
            # Constraint: Non-overlap with Sphere 1
            dist_to_s1 = np.linalg.norm(c_cand - c1)
            r_collision_s1 = dist_to_s1 - r1
            # Constraint: Non-overlap with adjacent corner sphere at (-x, x, z)
            dist_to_adj = 2 * x_coord
            r_collision_adj = dist_to_adj / 2.0
            
            r_potential = min(r_geom, r_collision_s1, r_collision_adj)
            
            if r_potential > best_r2:
                best_r2, best_c2_quadrant = r_potential, c_cand
    
    r2 = floor_to_step(best_r2, STEP)
    c2_base = np.round(best_c2_quadrant / STEP) * STEP
    
    # Add the 4 symmetric spheres
    for sx in [-1, 1]:
        for sy in [-1, 1]:
            c = np.array([c2_base[0] * sx, c2_base[1] * sy, c2_base[2]])
            spheres.append({'c': c, 'r': r2})

    # 4. Find Sphere 6 (On central axis, above Sphere 1)
    best_r6, best_c6 = 0, None
    z_start = c1[2] + r1 + STEP
    for z_coord in np.arange(z_start, H, STEP):
        c_cand = np.array([0, 0, z_coord])
        
        r_geom = max_radius_in_pyramid(c_cand)
        
        # Check collision with all 5 existing spheres
        r_collision = float('inf')
        for s in spheres:
            dist = np.linalg.norm(c_cand - s['c'])
            r_collision = min(r_collision, dist - s['r'])
        
        r_potential = min(r_geom, r_collision)
        
        if r_potential > best_r6:
            best_r6, best_c6 = r_potential, c_cand
            
    r6 = floor_to_step(best_r6, STEP)
    c6 = np.round(best_c6 / STEP) * STEP
    spheres.append({'c': c6, 'r': r6})

    # 5. Print results and find Min/Max radii
    print("Found optimal locations for 6 scans:\n")
    all_radii = []
    for i, s in enumerate(spheres):
        c, r = s['c'], s['r']
        # Final check on radius constraints
        final_r = max(MIN_R, min(r, MAX_R))
        all_radii.append(final_r)
        print(f"Scan {i+1}:")
        print(f"  Center (x, y, z) = ({c[0]:.1f}m, {c[1]:.1f}m, {c[2]:.1f}m)")
        print(f"  Radius r         = {final_r:.1f}m")

    max_radius = max(all_radii)
    min_radius = min(all_radii)
    
    print(f"\nMaximum scanning radius (R): {max_radius:.1f}m")
    print(f"Minimum scanning radius (r): {min_radius:.1f}m")
    
    final_answer = f"{max_radius:.1f}:{min_radius:.1f}"
    print(f"\nFinal Answer Format (R:r) = {final_answer}")
    return final_answer

final_answer_string = solve_pyramid_scanning()
print(f"<<<{final_answer_string}>>>")