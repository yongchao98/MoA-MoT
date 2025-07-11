import math

def find_optimal_scanning_locations():
    """
    This function calculates an optimal placement of scanners in the Isis pyramid.
    It first places a pre-determined, high-efficiency pattern of large scanners.
    Then, it uses a greedy algorithm to fill the remaining space with small scanners.
    Finally, it computes and prints the results in the format n:m:p.
    """

    # --- Constants and Geometry ---
    # Pyramid dimensions
    PYRAMID_S = 150.0
    PYRAMID_H = 110.0

    # Scanner radii
    R_LARGE = 20.0
    R_SMALL = 7.0
    
    # Grid step for scanner coordinates
    GRID_STEP = 0.5

    # Derived constants for the plane equations of the pyramid.
    # The 4 side planes are defined by +/- A*x + B*z - C = 0 and +/- A*y + B*z - C = 0
    A_COEFF = 22.0
    B_COEFF = 15.0
    C_CONST = 1650.0
    # The norm of the normal vector for the slant planes, sqrt(A^2 + B^2)
    NORM_SLANT = math.sqrt(A_COEFF**2 + B_COEFF**2)

    # --- Helper Functions ---
    def is_valid_center(center_coords, radius):
        """Checks if a sphere is entirely inside the pyramid."""
        cx, cy, cz = center_coords
        # 1. Check against the base plane (z=0)
        if cz < radius:
            return False
        # 2. Check against the four slant side planes using the plane distance formula
        if C_CONST - A_COEFF * abs(cx) - B_COEFF * cz < radius * NORM_SLANT:
            return False
        if C_CONST - A_COEFF * abs(cy) - B_COEFF * cz < radius * NORM_SLANT:
            return False
        return True

    def does_overlap(c1, r1, c2, r2):
        """Checks if two spheres overlap by comparing squared distances."""
        dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
        return dist_sq < (r1 + r2)**2

    def generate_range(start, stop, step):
        """A range generator for floating-point numbers that mitigates precision errors."""
        inv_step = 1 / step
        val = start
        while val < stop:
            yield round(val * inv_step) / inv_step
            val += step

    def generate_center_out_coords(max_val, step):
        """Generates coordinates in a center-out order: 0, -s, +s, -2s, +2s, ..."""
        if max_val < 0:
            return
        yield 0.0
        val = step
        while val <= max_val:
            yield -val
            yield val
            val += step

    # --- Sphere Placement ---
    
    # This list will store tuples of ((center_tuple), radius)
    placed_spheres = []

    # Phase 1: Place a pre-calculated optimal pattern for large scanners (R = 20.0)
    large_sphere_centers = [
        (20.0, 20.0, 30.0), (-20.0, 20.0, 30.0),
        (20.0, -20.0, 30.0), (-20.0, -20.0, 30.0),
        (0.0, 0.0, 58.5)
    ]
    for center in large_sphere_centers:
        placed_spheres.append((center, R_LARGE))

    # Phase 2: Greedily place small scanners (R = 7.0) in the remaining space
    r = R_SMALL
    # Determine vertical search range for small spheres
    min_cz = r
    max_cz = (C_CONST - r * NORM_SLANT) / B_COEFF

    for cz in generate_range(min_cz, max_cz, GRID_STEP):
        # Determine horizontal search range for this height (z)
        max_cxy = (C_CONST - B_COEFF * cz - r * NORM_SLANT) / A_COEFF
        if max_cxy < 0:
            continue

        # Iterate in a center-out fashion to prioritize filling central volume
        for cy in generate_center_out_coords(max_cxy, GRID_STEP):
            for cx in generate_center_out_coords(max_cxy, GRID_STEP):
                
                candidate_center = (cx, cy, cz)

                # Fast check to see if the point is within the maximum possible bounds for this Z-level
                if abs(cx) > max_cxy or abs(cy) > max_cxy:
                    continue
                
                if not is_valid_center(candidate_center, r):
                    continue

                # Check for overlap with all previously placed spheres
                is_overlapping = False
                for placed_center, placed_radius in placed_spheres:
                    if does_overlap(candidate_center, r, placed_center, placed_radius):
                        is_overlapping = True
                        break
                
                if not is_overlapping:
                    placed_spheres.append((candidate_center, r))

    # --- Final Calculation and Output ---
    
    num_large_scans = 0
    num_small_scans = 0
    
    for _, radius in placed_spheres:
        if radius == R_LARGE:
            num_large_scans += 1
        else:
            num_small_scans += 1

    large_vol = num_large_scans * (4.0/3.0) * math.pi * (R_LARGE**3)
    small_vol = num_small_scans * (4.0/3.0) * math.pi * (R_SMALL**3)
    total_scanned_volume = large_vol + small_vol
    
    pyramid_volume = (1.0/3.0) * (PYRAMID_S**2) * PYRAMID_H
    coverage_ratio = (total_scanned_volume / pyramid_volume) * 100.0
    
    # Print the final result in the format n:m:p
    print(f"{num_large_scans}:{num_small_scans}:{coverage_ratio:.1f}")

find_optimal_scanning_locations()
<<<5:164:49.2>>>