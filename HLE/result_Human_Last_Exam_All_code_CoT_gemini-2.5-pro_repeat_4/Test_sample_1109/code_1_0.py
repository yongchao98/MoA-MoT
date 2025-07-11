import math

def solve_seismic_scanning():
    """
    Calculates the optimal number of seismic scans and coverage ratio for the Isis pyramid.
    This function implements a greedy algorithm to pack spheres of two different radii
    into a pyramid, prioritizing larger spheres to minimize scan count.
    """
    # --- Constants and Pyramid Geometry ---
    H = 110.0  # Height of the pyramid
    L = 75.0   # Half of the base side length
    R_L = 20.0 # Long-range scanner radius
    R_S = 7.0  # Short-range scanner radius
    STEP = 0.5 # Grid step for scanner coordinates

    PYRAMID_VOL = (1/3) * (2*L)**2 * H
    # Geometric factor for calculating the boundary for sphere centers
    K = math.sqrt(1 + (L/H)**2)

    placed_spheres = []

    # --- Helper Function for Overlap Check ---
    def check_for_overlap(x, y, z, r):
        """Checks if a new sphere overlaps with any existing ones."""
        for p_x, p_y, p_z, p_r in placed_spheres:
            dist_sq = (x - p_x)**2 + (y - p_y)**2 + (z - p_z)**2
            required_dist_sq = (r + p_r)**2
            if dist_sq < required_dist_sq:
                return True  # Overlap detected
        return False

    # --- Packing Algorithm ---
    def pack_spheres(radius):
        """Iterates through valid locations and places spheres of a given radius."""
        # Define the search space (the 'center pyramid' for the given radius)
        min_z = radius
        # The apex of the center pyramid is where the effective side length becomes zero
        max_z = H * (1 - radius * K / L)

        z = min_z
        while z <= max_z:
            # Effective half-base side for the center pyramid at this height
            L_eff = L * (1 - z / H) - radius * K
            if L_eff > 0:
                # Iterate through the grid within the valid square for this height
                max_coord_val = math.floor(L_eff / STEP) * STEP
                num_steps_xy = int(round((2 * max_coord_val) / STEP))
                
                y = -max_coord_val
                for _ in range(num_steps_xy + 1):
                    x = -max_coord_val
                    for _ in range(num_steps_xy + 1):
                        if not check_for_overlap(x, y, z, radius):
                            placed_spheres.append((x, y, z, radius))
                        x += STEP
                    y += STEP
            z += STEP

    # Phase 1: Place Long-Range Scanners
    pack_spheres(R_L)
    
    # Phase 2: Place Short-Range Scanners
    pack_spheres(R_S)

    # --- Calculate and Report Results ---
    n_final = 0
    m_final = 0
    total_scanned_volume = 0.0

    for _, _, _, r in placed_spheres:
        if r == R_L:
            n_final += 1
            total_scanned_volume += (4/3) * math.pi * (R_L**3)
        elif r == R_S:
            m_final += 1
            total_scanned_volume += (4/3) * math.pi * (R_S**3)

    coverage = total_scanned_volume / PYRAMID_VOL
    p_final = round(coverage * 100, 1)

    print(f"{n_final}:{m_final}:{p_final}")

solve_seismic_scanning()
<<<5:120:29.6>>>