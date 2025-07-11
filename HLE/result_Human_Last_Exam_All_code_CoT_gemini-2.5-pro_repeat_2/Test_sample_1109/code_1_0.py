import math

def solve_seismic_scanning():
    """
    Calculates an optimal placement of seismic scanners in a pyramid
    using a greedy packing algorithm.
    """
    # 1. Define constants and parameters
    A = 150.0  # Pyramid base side length
    H = 110.0  # Pyramid height
    GRID_STEP = 0.5

    R_LONG = 20.0
    R_SHORT = 7.0

    # Derived constants for calculations
    PI = math.pi
    V_PYRAMID = (1.0/3.0) * A**2 * H
    SQRT_TERM = math.sqrt(4 * H**2 + A**2)
    AH = A * H

    # List to store placed scanners: (x, y, z, R)
    placed_scans = []

    # 2. Helper function to check for overlaps
    def is_not_overlapping(cx, cy, cz, r_new, existing_scans):
        for px, py, pz, pr_old in existing_scans:
            dist_sq = (cx - px)**2 + (cy - py)**2 + (cz - pz)**2
            min_dist = r_new + pr_old
            # Use a small tolerance for floating point comparisons
            if dist_sq < min_dist**2 - 1e-9:
                return False
        return True

    # 3. Implement the greedy packing algorithm
    scanner_radii = [R_LONG, R_SHORT]

    for r_scan in scanner_radii:
        # Define the search space for the center of the sphere (the "inner pyramid")
        bound = AH - r_scan * SQRT_TERM
        
        # Determine the valid range for the z-coordinate of the center
        cz_min = r_scan
        cz_max = bound / A # Max z if cx and cy are 0
        
        # Use integer-based loops for numerical stability
        z_k_min = math.ceil(cz_min / GRID_STEP)
        z_k_max = math.floor(cz_max / GRID_STEP)

        for k in range(int(z_k_min), int(z_k_max) + 1):
            cz = k * GRID_STEP
            
            # Determine the valid x/y range for the current z
            xy_max_at_z = (bound - A * cz) / (2 * H)
            if xy_max_at_z < 0:
                continue
            
            xy_j_max = math.floor(xy_max_at_z / GRID_STEP)

            for j in range(-int(xy_j_max), int(xy_j_max) + 1):
                cy = j * GRID_STEP
                for i in range(-int(xy_j_max), int(xy_j_max) + 1):
                    cx = i * GRID_STEP
                    
                    # Check for overlaps with already placed scanners
                    if is_not_overlapping(cx, cy, cz, r_scan, placed_scans):
                        placed_scans.append((cx, cy, cz, r_scan))

    # 4. Calculate final results
    n_long = 0
    m_short = 0
    v_scanned = 0.0

    for _, _, _, r in placed_scans:
        if r == R_LONG:
            n_long += 1
        elif r == R_SHORT:
            m_short += 1
    
    v_scanned_long = n_long * (4.0/3.0) * PI * R_LONG**3
    v_scanned_short = m_short * (4.0/3.0) * PI * R_SHORT**3
    v_scanned = v_scanned_long + v_scanned_short

    coverage_ratio = (v_scanned / V_PYRAMID) * 100

    # Output the final answer in the required format n:m:p
    print(f"{n_long}:{m_short}:{coverage_ratio:.1f}")

if __name__ == '__main__':
    solve_seismic_scanning()
<<<13:210:35.3>>>