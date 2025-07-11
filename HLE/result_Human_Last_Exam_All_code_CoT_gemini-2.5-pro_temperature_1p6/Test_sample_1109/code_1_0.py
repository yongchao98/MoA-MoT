import math

def solve_seismic_scanning():
    """
    Calculates the optimal placement of seismic scanners in the Isis pyramid.
    This script uses a greedy algorithm to pack spheres (scans) into the pyramid,
    prioritizing larger scans and filling from the bottom up.
    """
    # 1. Define Constants and Geometry
    
    # Pyramid dimensions
    H = 110.0  # height in meters
    A = 75.0   # half of the base side length in meters
    
    # Scanner parameters
    R_LONG = 20.0  # long-range radius
    R_SHORT = 7.0  # short-range radius
    STEP = 0.5     # grid step for scanner coordinates
    
    # Pre-calculated geometric and volumetric values for efficiency
    V_PYRAMID = (1/3) * (2*A)**2 * H
    V_LONG = (4/3) * math.pi * R_LONG**3
    V_SHORT = (4/3) * math.pi * R_SHORT**3
    
    # This value is used in the sphere containment check. It's derived from the
    # formula for the distance from a point to a plane.
    SQRT_H2_A2 = math.sqrt(H**2 + A**2)
    
    # The condition for a sphere to be inside is: H*max(|cx|,|cy|) + A*cz <= A*H - R*sqrt(H^2+A^2)
    # We pre-calculate the right side of the inequality for each radius.
    LIMIT_LONG = A * H - R_LONG * SQRT_H2_A2
    LIMIT_SHORT = A * H - R_SHORT * SQRT_H2_A2

    # 2. Helper Functions
    
    def is_center_valid(p, R, limit_val):
        """Checks if a sphere with center p and radius R is fully inside the pyramid."""
        cx, cy, cz = p
        # Sphere must be above the base plane (z=0) by at least its radius.
        if cz < R:
            return False
        # Check if the sphere is inside the four side planes.
        if H * max(abs(cx), abs(cy)) + A * cz > limit_val:
            return False
        return True

    def check_overlap(p1, r1, placed_spheres):
        """Checks if a new sphere (p1, r1) overlaps with any already placed spheres."""
        for p2, r2 in placed_spheres:
            # Using squared distances is faster as it avoids square roots.
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            if dist_sq < (r1 + r2)**2:
                return True
        return False

    # 3. Greedy Packing Algorithm
    
    placed_spheres = []
    
    # Generate the grid coordinates, sorted in a center-out fashion
    max_abs_coord = int(A / STEP)
    coord_vals = [i * STEP for i in range(-max_abs_coord, max_abs_coord + 1)]
    xy_coords = sorted(coord_vals, key=lambda v: (abs(v), v))
    
    # Phase 1: Place Long-Range Scanners
    z_coords_long = (i * STEP for i in range(int(R_LONG / STEP), int(H / STEP)))
    for cz in z_coords_long:
        # Optimization: determine max valid x/y for this height to reduce search space
        max_c = (LIMIT_LONG - A * cz) / H
        if max_c < 0: continue
        for cx in xy_coords:
            if abs(cx) > max_c: break
            for cy in xy_coords:
                if abs(cy) > max_c: break
                
                center = (cx, cy, cz)
                if not check_overlap(center, R_LONG, placed_spheres):
                    placed_spheres.append((center, R_LONG))

    n_long = len(placed_spheres)
    
    # Phase 2: Place Short-Range Scanners
    z_coords_short = (i * STEP for i in range(int(R_SHORT / STEP), int(H / STEP)))
    for cz in z_coords_short:
        max_c = (LIMIT_SHORT - A * cz) / H
        if max_c < 0: continue
        for cx in xy_coords:
            if abs(cx) > max_c: break
            for cy in xy_coords:
                if abs(cy) > max_c: break
                
                center = (cx, cy, cz)
                # We still must do a full validity check because the loop optimization
                # only checks one dimension of max(|cx|, |cy|).
                if is_center_valid(center, R_SHORT, LIMIT_SHORT):
                    if not check_overlap(center, R_SHORT, placed_spheres):
                        placed_spheres.append((center, R_SHORT))
                        
    n_short = len(placed_spheres) - n_long

    # 4. Calculate and Print Final Metrics
    total_scanned_volume = n_long * V_LONG + n_short * V_SHORT
    coverage_p = (total_scanned_volume / V_PYRAMID) * 100.0

    print("--- Seismic Scanning Result ---")
    print(f"An optimal placement strategy resulted in:")
    print(f"  n (long range scans): {n_long}")
    print(f"  m (short range scans): {n_short}")
    print("\nCalculation of Coverage Ratio (p):")
    print(f"  Volume of one long range scan: {V_LONG:.1f} m^3")
    print(f"  Volume of one short range scan: {V_SHORT:.1f} m^3")
    print(f"  Total Pyramid Volume: {V_PYRAMID:.1f} m^3")
    print(f"\n  Equation:")
    print(f"  (({n_long} * {V_LONG:.1f}) + ({n_short} * {V_SHORT:.1f})) / {V_PYRAMID:.1f} = {coverage_p:.1f}%")
    
    print("\nFinal Answer (n:m:p):")
    print(f"{n_long}:{n_short}:{coverage_p:.1f}")


solve_seismic_scanning()