import math
import sys

def solve_pyramid_packing():
    """
    This script finds the optimal placement of 6 spheres inside a pyramid
    to maximize their total volume using a grid search on a symmetric configuration.
    
    The configuration consists of:
    - 1 sphere on the central z-axis (bottom).
    - 4 spheres in a square pattern on a horizontal plane.
    - 1 sphere on the central z-axis (top).
    
    The script iterates through possible positions, calculates the maximum valid radii 
    for each position, and finds the arrangement with the largest total volume.
    """
    
    # --- Problem Constants ---
    H = 110.0  # Pyramid height
    A = 150.0  # Pyramid base side length
    MIN_R = 10.0
    MAX_R = 50.0
    STEP = 0.5 # Discretization step for coordinates and radii

    # Derived geometric constants for the pyramid's slanted faces
    K = H / (A / 2.0)
    NORM_FACTOR = math.sqrt(K**2 + 1)

    # --- Search Initialization ---
    max_volume = 0
    best_config = {}
    
    # --- Search Ranges ---
    # These ranges are chosen to be broad enough to find a good solution 
    # while keeping the execution time manageable.
    cz_below_range = [i * STEP for i in range(int(10/STEP), int(25/STEP))]
    cz_plane_range = [i * STEP for i in range(int(20/STEP), int(50/STEP))]
    d_range =        [i * STEP for i in range(int(15/STEP), int(45/STEP))]

    # --- Search Execution ---
    print("Searching for the optimal scanner configuration...")
    print("This process may take a few minutes depending on your computer's speed.")
    
    total_iterations = len(cz_below_range) * len(cz_plane_range) * len(d_range)
    iter_count = 0
    
    # Main search loop iterates over possible positions for the bottom and plane spheres.
    for cz_below in cz_below_range:
        for cz_plane in cz_plane_range:
            for d in d_range:
                # Progress indicator
                iter_count += 1
                if iter_count % 10000 == 0:
                    progress = (iter_count / total_iterations) * 100
                    sys.stdout.write(f"\rProgress: {progress:.1f}%")
                    sys.stdout.flush()

                # --- Sphere Below ---
                # Calculate max possible radius based on its position
                r_below_max = min(cz_below, (H - cz_below) / NORM_FACTOR)
                r_below = math.floor(r_below_max * (1/STEP)) / (1/STEP)
                if not (MIN_R <= r_below <= MAX_R): continue
                
                # --- Plane Spheres ---
                # Calculate max possible radius based on their position
                r_plane_max = min(d, cz_plane, (H - K * d - cz_plane) / NORM_FACTOR)
                r_plane = math.floor(r_plane_max * (1/STEP)) / (1/STEP)
                if not (MIN_R <= r_plane <= MAX_R): continue
                
                # Check for overlap between the plane spheres and the bottom sphere
                dist_sq_pb = 2 * d**2 + (cz_plane - cz_below)**2
                if dist_sq_pb < (r_plane + r_below)**2: continue

                # --- Sphere Above ---
                # Now, find the best possible top sphere for this base configuration
                # We search for its z-coordinate, starting just above the plane spheres
                cz_above_start = cz_plane + r_plane + MIN_R 
                
                # Iterate upwards to find an optimal position for the top sphere
                for cz_above in [i * STEP for i in range(int(cz_above_start / STEP), int(100 / STEP))]:
                    r_above_max = min(cz_above, (H - cz_above) / NORM_FACTOR)
                    r_above = math.floor(r_above_max * (1/STEP)) / (1/STEP)
                    if not (MIN_R <= r_above <= MAX_R): continue
                    
                    # Check for overlaps with the top sphere
                    if (cz_above - cz_below) < (r_above + r_below): continue
                    
                    dist_sq_pa = 2 * d**2 + (cz_above - cz_plane)**2
                    if dist_sq_pa < (r_plane + r_above)**2: continue

                    # If we are here, we have a valid 6-sphere configuration
                    current_volume = 4 * r_plane**3 + r_below**3 + r_above**3
                    if current_volume > max_volume:
                        max_volume = current_volume
                        best_config = {
                            'volume': current_volume,
                            'below': {'c': (0, 0, cz_below), 'r': r_below},
                            'plane': {'d': d, 'r': r_plane, 'cz': cz_plane},
                            'above': {'c': (0, 0, cz_above), 'r': r_above}
                        }
                    # Since r_above will only decrease as cz_above increases from here,
                    # we can stop searching for a better top sphere for this base.
                    break 

    # --- Print Results ---
    sys.stdout.write("\rSearch complete. Finalizing results...\n")
    
    if not best_config:
        print("No valid configuration found with the given search parameters.")
    else:
        # Unpack the best found configuration
        r_b = best_config['below']['r']
        cz_b = best_config['below']['c'][2]
        r_p = best_config['plane']['r']
        d_p = best_config['plane']['d']
        cz_p = best_config['plane']['cz']
        r_a = best_config['above']['r']
        cz_a = best_config['above']['c'][2]

        all_radii = [r_b, r_p, r_a]
        max_r = max(all_radii)
        min_r = min(all_radii)

        print("\n--- Optimal Configuration Found ---")
        print(f"Maximized Volume Factor (sum of r^3): {best_config['volume']:.1f}")
        print("\nDetails of the 6 Scans:")
        print("  - 1 Bottom Sphere: Radius={:.1f}m, Center=(0.0, 0.0, {:.1f})m".format(r_b, cz_b))
        print("  - 4 Plane Spheres: Radius={:.1f}m, Centers=(+/i-{:.1f}, +/-{:.1f}, {:.1f})m".format(r_p, d_p, d_p, cz_p))
        print("  - 1 Top Sphere:    Radius={:.1f}m, Center=(0.0, 0.0, {:.1f})m".format(r_a, cz_a))

        print("\nMaximum and minimum scanning radii (R:r):")
        # Final answer format
        print(f"<<<{max_r:.1f}:{min_r:.1f}>>>")

if __name__ == '__main__':
    solve_pyramid_packing()