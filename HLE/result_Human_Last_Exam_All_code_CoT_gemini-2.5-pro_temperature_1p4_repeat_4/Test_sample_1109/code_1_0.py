import math

def solve_seismic_scanning():
    """
    This function calculates the optimal placement of seismic scanners in the Isis pyramid
    based on a greedy algorithm, prioritizing larger scanners first to maximize coverage
    while minimizing the number of scans.
    """
    # Step 1: Define constants and parameters for the pyramid and scanners.
    BASE_SIDE = 150.0
    HEIGHT = 110.0
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5
    
    # Pre-calculate geometric constants for efficiency.
    # ALPHA is related to the slope of the pyramid's side faces.
    ALPHA = (BASE_SIDE / 2.0) / HEIGHT
    # BETA is a geometric factor used to define the valid region for sphere centers.
    BETA = math.sqrt(1 + ALPHA**2)
    
    # Calculate the total volume of the pyramid.
    V_PYRAMID = (1.0/3.0) * (BASE_SIDE**2) * HEIGHT

    def sphere_volume(R):
        """Calculates the volume of a sphere with radius R."""
        return (4.0/3.0) * math.pi * R**3

    def get_placements(R, existing_spheres):
        """
        Finds valid, non-overlapping placements for spheres of a given radius.
        It uses a greedy approach, iterating from the bottom-up and center-out.
        To handle the pyramid's symmetry efficiently, it calculates placements
        in one quadrant and then mirrors the results.
        """
        new_placements = []
        
        # Calculate the valid z-range for the center of a sphere of radius R.
        cz_min = R
        cz_max = ((BASE_SIDE / 2.0) - R * BETA) / ALPHA
        
        # Create a list of z-coordinates to check, respecting the grid step.
        z_coords = [i * GRID_STEP for i in range(math.ceil(cz_min / GRID_STEP), int(cz_max / GRID_STEP) + 1)]

        for cz in z_coords:
            # For each height `cz`, calculate the maximum horizontal distance `W_cz` from the
            # central axis where a sphere's center can be placed.
            W_cz = (BASE_SIDE / 2.0) - ALPHA * cz - R * BETA
            if W_cz < 0:
                continue
            
            # Create a list of x/y coordinates to check in the first quadrant (x>=0, y>=0).
            c_limit = math.floor(W_cz / GRID_STEP) * GRID_STEP
            c_coords_q1 = [i * GRID_STEP for i in range(int(c_limit / GRID_STEP) + 1)]

            for cx in c_coords_q1:
                for cy in c_coords_q1:
                    # Use a set to generate symmetric points, automatically handling duplicates on axes.
                    symmetric_points = set()
                    symmetric_points.add((cx, cy, cz))
                    symmetric_points.add((-cx, cy, cz))
                    symmetric_points.add((cx, -cy, cz))
                    symmetric_points.add((-cx, -cy, cz))
                    
                    # Obstacles include spheres from previous passes and new ones from this pass.
                    all_obstacles = existing_spheres + new_placements

                    # Sort for deterministic behavior
                    for p_cand in sorted(list(symmetric_points)):
                        is_valid = True
                        for obs_center, obs_R in all_obstacles:
                            # Use squared distances to avoid costly sqrt operations.
                            dist_sq = (p_cand[0] - obs_center[0])**2 + \
                                      (p_cand[1] - obs_center[1])**2 + \
                                      (p_cand[2] - obs_center[2])**2
                            min_dist_sq = (R + obs_R)**2
                            if dist_sq < min_dist_sq:
                                is_valid = False
                                break
                        
                        if is_valid:
                            new_placements.append((p_cand, R))
                            # Update the list of obstacles for checking the remaining symmetric points.
                            all_obstacles.append((p_cand, R))
                            
        return new_placements

    print("Calculating optimal scanner locations...\n")
    
    # List to store all placed spheres as ((center_x, y, z), radius)
    placed_spheres = []

    # Place long-range scanners.
    long_placements = get_placements(R_LONG, placed_spheres)
    placed_spheres.extend(long_placements)
    n = len(long_placements)

    # Place short-range scanners in the remaining space.
    short_placements = get_placements(R_SHORT, placed_spheres)
    m = len(short_placements)

    # Calculate the final results.
    v_long = sphere_volume(R_LONG)
    v_short = sphere_volume(R_SHORT)
    total_scanned_volume = n * v_long + m * v_short
    coverage_ratio = (total_scanned_volume / V_PYRAMID) * 100

    # Output the results, including the equation as requested.
    print("--- Calculation Details ---")
    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print(f"Volume of one long-range scan: {v_long:.2f} m^3")
    print(f"Volume of one short-range scan: {v_short:.2f} m^3")
    print(f"Total volume of the pyramid: {V_PYRAMID:.2f} m^3\n")
    
    print("Coverage Ratio Equation:")
    print(f"p = (n * Volume_long + m * Volume_short) / Volume_pyramid")
    print(f"p = ({n} * {v_long:.2f} + {m} * {v_short:.2f}) / {V_PYRAMID:.2f}")
    print(f"p = ({n * v_long:.2f} + {m * v_short:.2f}) / {V_PYRAMID:.2f}")
    print(f"p = {total_scanned_volume:.2f} / {V_PYRAMID:.2f} = {total_scanned_volume/V_PYRAMID:.5f}")
    print(f"Coverage Ratio (p) = {coverage_ratio:.1f}%\n")

    final_answer = f"{n}:{m}:{coverage_ratio:.1f}"
    print("--- Final Answer ---")
    print(f"<<<{final_answer}>>>")

solve_seismic_scanning()