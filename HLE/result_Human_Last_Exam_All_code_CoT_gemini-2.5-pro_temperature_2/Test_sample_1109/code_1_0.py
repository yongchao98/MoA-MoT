import math

def solve_seismic_scanning():
    """
    Finds an optimal placement of seismic scanners in a pyramid using a greedy algorithm.
    """
    # --- Constants and Setup ---
    PYRAMID_H = 110.0
    PYRAMID_A = 150.0
    R_LONG = 20.0
    R_SHORT = 7.0
    SCAN_STEP = 0.5

    V_PYRAMID = (1/3) * (PYRAMID_A**2) * PYRAMID_H
    V_LONG = (4/3) * math.pi * (R_LONG**3)
    V_SHORT = (4/3) * math.pi * (R_SHORT**3)
    
    # Pre-calculated constant for distance formula to slanted pyramid faces
    SQRT_709 = math.sqrt(22**2 + 15**2)

    placed_spheres = []

    # --- Helper Functions ---
    def check_overlap(new_center, new_radius, existing_spheres):
        """Checks if a new sphere overlaps with any existing ones."""
        ncx, ncy, ncz = new_center
        for s in existing_spheres:
            scx, scy, scz = s['center']
            # Using squared distances to avoid costly square root operations
            dist_sq = (ncx - scx)**2 + (ncy - scy)**2 + (ncz - scz)**2
            min_dist = new_radius + s['radius']
            if dist_sq < min_dist**2:
                return True
        return False

    def generate_ordered_coords(max_val, step):
        """Generates coordinates spreading from the center outwards."""
        coords = [0.0]
        num_steps = int(max_val / step)
        for i in range(1, num_steps + 1):
            val = i * step
            coords.append(-val)
            coords.append(val)
        return coords

    # --- Packing Loop Function ---
    def pack_one_type(radius, existing_spheres):
        """
        Greedily places scanners of a given radius.
        Returns a list of newly placed spheres.
        """
        new_spheres = []
        # Define the geometric boundaries for valid centers
        z_min = radius
        z_max = 110.0 - (radius * SQRT_709 / 15.0)

        # Iterate from the top of the valid region downwards for Z
        z_coords_int = range(int(z_min * 2), int(z_max * 2) + 1)
        for zc_int in sorted(list(z_coords_int), reverse=True):
            zc = zc_int / 2.0
            
            # For each Z, determine the max horizontal extent for a center
            w_max = (1650 - radius * SQRT_709 - 15 * zc) / 22.0
            if w_max < 0:
                continue

            # Iterate X and Y from the center outwards
            for yc in generate_ordered_coords(w_max, SCAN_STEP):
                for xc in generate_ordered_coords(w_max, SCAN_STEP):
                    # Check if point is within the tighter horizontal boundary
                    if abs(xc) > w_max or abs(yc) > w_max:
                        continue
                    
                    candidate_center = (xc, yc, zc)
                    
                    # Check for overlap with all spheres placed so far
                    if not check_overlap(candidate_center, radius, existing_spheres + new_spheres):
                        new_spheres.append({'center': candidate_center, 'radius': radius})
        return new_spheres

    # --- Main Execution ---
    # Phase 1: Place long-range scanners
    long_range_spheres = pack_one_type(R_LONG, placed_spheres)
    placed_spheres.extend(long_range_spheres)

    # Phase 2: Fill gaps with short-range scanners
    short_range_spheres = pack_one_type(R_SHORT, placed_spheres)
    placed_spheres.extend(short_range_spheres)
    
    # --- Final Calculation and Output ---
    n = len(long_range_spheres)
    m = len(short_range_spheres)

    total_scanned_volume = n * V_LONG + m * V_SHORT
    p = (total_scanned_volume / V_PYRAMID) * 100
    
    # Print the result in the specified format "n:m:p"
    print(f"{n}:{m}:{p:.1f}")

solve_seismic_scanning()
<<<8:188:65.2>>>