import math
import numpy as np

def solve_seismic_scanning():
    """
    Finds an optimal set of scanning locations inside a pyramid.
    """
    # 1. Define Geometry and Constants
    L = 150.0  # Base side in meters
    H = 110.0  # Height in meters
    R_L = 20.0 # Long range radius
    R_S = 7.0  # Short range radius
    GRID_STEP = 0.5

    # Pre-calculated values for performance
    V_PYRAMID = (1/3) * L**2 * H
    V_SPHERE_L = (4/3) * math.pi * R_L**3
    V_SPHERE_S = (4/3) * math.pi * R_S**3
    SIDE_PLANE_DIST_FACTOR = math.sqrt(4 * H**2 + L**2)
    L_H = L * H

    def is_valid_center(c, r):
        """Checks if a sphere (center c, radius r) is inside the pyramid."""
        cx, cy, cz = c
        if cz < r:
            return False
        # Condition derived from plane distance formula
        # 2H*|c| + L*cz - LH + r*sqrt(4H^2+L^2) <= 0
        if 2 * H * abs(cx) + L * cz - L_H + r * SIDE_PLANE_DIST_FACTOR > 1e-9:
           return False
        if 2 * H * abs(cy) + L * cz - L_H + r * SIDE_PLANE_DIST_FACTOR > 1e-9:
           return False
        return True

    def generate_candidates(r):
        """Generates all valid candidate centers for a given radius on the grid."""
        candidates = []
        z_coords = np.arange(r, H, GRID_STEP)
        for cz in z_coords:
            max_pyr_coord = (L / 2.0) * (1.0 - cz / H)
            if max_pyr_coord <= 0:
                continue
            
            # Iterate a bounding box for the pyramid's cross-section
            xy_limit = int(max_pyr_coord / GRID_STEP)
            for ix in range(-xy_limit, xy_limit + 1):
                cx = ix * GRID_STEP
                for iy in range(-xy_limit, xy_limit + 1):
                    cy = iy * GRID_STEP
                    if is_valid_center((cx, cy, cz), r):
                        candidates.append((cx, cy, cz))
        return candidates
        
    def dist_sq(p1, p2):
        """Calculates the squared distance between two points."""
        return (p1[0]-p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]-p2[2])**2

    # --- Main Logic ---
    placed_scans = []
    
    # 2. & 3. Place Long-Range Scanners
    print("Placing long-range scanners (r=20m)...")
    candidates_L = generate_candidates(R_L)
    candidates_L.sort(key=lambda p: (p[2], abs(p[0]), abs(p[1])))
    
    min_dist_LL_sq = (2 * R_L)**2
    while candidates_L:
        p_new = candidates_L[0]
        placed_scans.append({'center': p_new, 'radius': R_L})
        candidates_L = [p for p in candidates_L if dist_sq(p, p_new) >= min_dist_LL_sq]
    
    n = len(placed_scans)
    print(f"Placed {n} long-range scanners.")

    # 4. Place Short-Range Scanners
    print("Placing short-range scanners (r=7m)...")
    candidates_S = generate_candidates(R_S)
    
    # Filter against already placed long-range scans
    min_dist_LS_sq = (R_L + R_S)**2
    for scan_L in placed_scans:
        center_L = scan_L['center']
        candidates_S = [p for p in candidates_S if dist_sq(p, center_L) >= min_dist_LS_sq]
    
    candidates_S.sort(key=lambda p: (p[2], abs(p[0]), abs(p[1])))
    
    min_dist_SS_sq = (2 * R_S)**2
    num_L_scans = len(placed_scans)
    while candidates_S:
        p_new = candidates_S[0]
        placed_scans.append({'center': p_new, 'radius': R_S})
        candidates_S = [p for p in candidates_S if dist_sq(p, p_new) >= min_dist_SS_sq]

    # 5. Calculate and Print Final Result
    m = len(placed_scans) - num_L_scans
    print(f"Placed {m} short-range scanners.")

    v_scanned = n * V_SPHERE_L + m * V_SPHERE_S
    coverage_ratio = (v_scanned / V_PYRAMID) * 100
    
    print("\n--- Final Result ---")
    # We are asked to output each number in the final equation.
    # The variables are n, m, p.
    # The output format is n:m:p
    print(f"{n}:{m}:{coverage_ratio:.1f}")
    
solve_seismic_scanning()