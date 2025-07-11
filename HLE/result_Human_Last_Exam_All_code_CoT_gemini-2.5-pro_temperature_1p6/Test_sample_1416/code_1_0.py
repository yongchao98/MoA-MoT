import math
import random

def solve_pyramid_scanning():
    """
    Finds optimal locations for 6 seismic scans in a pyramid using a greedy search algorithm.
    """
    # Define pyramid and scanner parameters
    PYRAMID_H = 110.0
    PYRAMID_S = 150.0
    PYRAMID_B = PYRAMID_S / 2.0
    R_MIN, R_MAX = 10.0, 50.0
    STEP = 0.5
    N_SCANS = 6
    
    # Use a fixed seed for reproducible results
    random.seed(3)
    
    # Number of random candidate centers to evaluate for each scan
    N_CANDIDATES = 200000 
    
    # Pre-calculated constant for distance to slanted planes
    K = math.sqrt(PYRAMID_H**2 + PYRAMID_B**2)

    def round_to_step(value, step):
        """Rounds a value down to the nearest multiple of the step."""
        return math.floor(value / step) * step

    def get_max_radius_in_pyramid(x, y, z):
        """Calculates the max radius of a sphere at (x, y, z) that fits in the pyramid."""
        x_abs, y_abs = abs(x), abs(y)
        
        # A point is outside the pyramid's cross-section at height z
        if x_abs > PYRAMID_B * (1 - z / PYRAMID_H) or y_abs > PYRAMID_B * (1 - z / PYRAMID_H):
            return 0

        # Distance to base plane (z=0)
        dist_base = z
        
        # Distance to the four slanted side planes
        dist_side_x = (PYRAMID_H * PYRAMID_B - PYRAMID_H * x_abs - PYRAMID_B * z) / K
        dist_side_y = (PYRAMID_H * PYRAMID_B - PYRAMID_H * y_abs - PYRAMID_B * z) / K

        # Point must be inside, so distances should be non-negative
        if dist_side_x < 0 or dist_side_y < 0:
            return 0

        return min(dist_base, dist_side_x, dist_side_y)

    placed_spheres = []

    for i in range(N_SCANS):
        best_sphere_candidate = None
        max_r_cubed = -1

        for _ in range(N_CANDIDATES):
            # Generate a random candidate center point within the pyramid
            # Prioritize lower sections where larger spheres can fit
            z_c = random.uniform(R_MIN, PYRAMID_H - 0.01)
            max_coord_at_z = PYRAMID_B * (1 - z_c / PYRAMID_H)
            
            if max_coord_at_z < R_MIN:
                continue

            x_c = random.uniform(-max_coord_at_z, max_coord_at_z)
            y_c = random.uniform(-max_coord_at_z, max_coord_at_z)

            # Snap center coordinates to the required grid
            x_c = round_to_step(x_c, STEP)
            y_c = round_to_step(y_c, STEP)
            z_c = round_to_step(z_c, STEP)

            # --- Calculate max valid radius for this candidate center ---

            # 1. Start with the constraint from the pyramid walls
            r_potential = get_max_radius_in_pyramid(x_c, y_c, z_c)

            # 2. Reduce radius based on non-overlap with already placed spheres
            for sx, sy, sz, sr in placed_spheres:
                dist = math.sqrt((x_c - sx)**2 + (y_c - sy)**2 + (z_c - sz)**2)
                r_potential = min(r_potential, dist - sr)

            # 3. Apply scanner hardware constraints and snap to grid
            if r_potential < R_MIN:
                continue
            
            r_final = round_to_step(min(r_potential, R_MAX), STEP)

            if r_final < R_MIN:
                continue

            # Check if this candidate is the best found for the current scan
            if r_final**3 > max_r_cubed:
                max_r_cubed = r_final**3
                best_sphere_candidate = (x_c, y_c, z_c, r_final)

        if best_sphere_candidate:
            placed_spheres.append(best_sphere_candidate)
        else:
            # This indicates a failure to place a sphere, might happen if space is too tight.
            # We can simply report that we couldn't place all 6.
            pass

    # --- Print Results ---
    print("Optimal scanning locations and radii found for N=6:\n")
    all_radii = []
    
    if not placed_spheres:
        print("No valid scanning locations could be found.")
        return

    for i, (x, y, z, r) in enumerate(placed_spheres):
        all_radii.append(r)
        print(f"Scan {i+1}: Center=({x:.1f}, {y:.1f}, {z:.1f}) m, Radius={r:.1f} m")

    if len(placed_spheres) < N_SCANS:
        print(f"\nWarning: Only able to place {len(placed_spheres)} of {N_SCANS} scans.")

    max_r = max(all_radii)
    min_r = min(all_radii)
    
    print(f"\nFinal Result:")
    print(f"R:{max_r:.1f}, r:{min_r:.1f}")

solve_pyramid_scanning()