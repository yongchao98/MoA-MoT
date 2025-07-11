import numpy as np
import math

def solve():
    """
    Finds the optimal placement of 6 spheres inside a pyramid using a greedy algorithm.
    """
    # 1. Define Pyramid Geometry and Scanner Constraints
    PYRAMID_H = 110.0
    PYRAMID_A = 150.0
    MIN_R, MAX_R = 10.0, 50.0
    STEP = 0.5  # Coordinates and radii must be multiples of this
    N_SPHERES = 6
    
    # Pre-calculate for plane distance formula: 110x + 75z - 8250 = 0
    # The normal vector part is sqrt(110^2 + 75^2)
    PLANE_DENOMINATOR = math.sqrt(110**2 + 75**2)

    def get_max_radius_in_pyramid(c):
        """Calculates max radius for a sphere at center c to fit inside the pyramid."""
        cx, cy, cz = c
        
        # Check if center is outside the pyramid's horizontal slice
        max_width_at_cz = PYRAMID_A * (1 - cz / PYRAMID_H)
        if abs(cx) * 2 > max_width_at_cz or abs(cy) * 2 > max_width_at_cz or cz < 0 or cz > PYRAMID_H:
            return 0.0

        # Distance to base plane (z=0)
        r_base = cz
        
        # Distance to the four side planes
        # Equation for side plane in x>0 octant is 110*x + 75*z - 8250 = 0
        # Distance = |110*cx + 75*cz - 8250| / sqrt(110^2 + 75^2)
        # Since the center is inside, the term is negative, so we use -(...)
        r_sidex = (8250 - 110 * abs(cx) - 75 * cz) / PLANE_DENOMINATOR
        r_sidey = (8250 - 110 * abs(cy) - 75 * cz) / PLANE_DENOMINATOR
        
        return max(0, min(r_base, r_sidex, r_sidey))

    spheres = []
    
    # 2. Iteratively place N_SPHERES using a greedy approach
    for i in range(N_SPHERES):
        best_candidate = None
        max_r_cubed = -1.0
        
        # Define a search grid for the center of the next sphere
        # A coarser grid makes the search faster. 2.0m is a good compromise.
        search_step = 2.0 
        z_min_bound = MIN_R
        
        x_range = np.arange(-PYRAMID_A / 2 + search_step, PYRAMID_A / 2, search_step)
        y_range = np.arange(-PYRAMID_A / 2 + search_step, PYRAMID_A / 2, search_step)
        z_range = np.arange(z_min_bound, PYRAMID_H - search_step, search_step)
        
        for cz in z_range:
            for cy in y_range:
                for cx in x_range:
                    center_candidate = (cx, cy, cz)
                    
                    # Max radius allowed by pyramid walls
                    r_potential = get_max_radius_in_pyramid(center_candidate)
                    if r_potential < MIN_R:
                        continue
                        
                    # Max radius allowed by non-overlap with existing spheres
                    r_nonoverlap = float('inf')
                    for placed_sphere in spheres:
                        dist_to_placed = np.linalg.norm(np.array(center_candidate) - np.array(placed_sphere['center']))
                        r_nonoverlap = min(r_nonoverlap, dist_to_placed - placed_sphere['radius'])
                    
                    # Combine constraints
                    r_final = min(r_potential, r_nonoverlap)
                    
                    # Clamp to scanner's effective range
                    r_final = min(r_final, MAX_R)
                    
                    if r_final < MIN_R:
                        continue
                    
                    # Check if this sphere is better than the best one found so far
                    if r_final**3 > max_r_cubed:
                        max_r_cubed = r_final**3
                        best_candidate = {'center': center_candidate, 'radius': r_final}
        
        if best_candidate:
            # Refine coordinates and radius to be multiples of STEP
            best_r = math.floor(best_candidate['radius'] / STEP) * STEP
            
            # Ensure radius is still valid after rounding
            if best_r < MIN_R:
                continue

            best_c = tuple(round(coord / STEP) * STEP for coord in best_candidate['center'])
            
            # Add the best sphere found in this iteration to our list
            spheres.append({'center': best_c, 'radius': best_r})
        else:
            # If no more spheres can be placed, stop.
            print(f"Could only place {len(spheres)} spheres.")
            break

    # 3. Output the results
    print("Found optimal locations for 6 spheres:")
    all_radii = []
    if not spheres:
        print("Could not place any spheres meeting the criteria.")
        return

    for i, s in enumerate(spheres):
        c = s['center']
        r = s['radius']
        all_radii.append(r)
        print(f"Sphere {i+1}: Center=({c[0]:.1f}, {c[1]:.1f}, {c[2]:.1f}), Radius={r:.1f}")

    max_radius = max(all_radii)
    min_radius = min(all_radii)

    print(f"\nMax Radius (R): {max_radius:.1f}")
    print(f"Min Radius (r): {min_radius:.1f}")
    
    # Final answer format
    final_answer = f"{max_radius:.1f}:{min_radius:.1f}"
    print(f"\nFinal Answer (R:r): {final_answer}")
    print(f"<<<{final_answer}>>>")


solve()