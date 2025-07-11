import math
import numpy as np

def solve():
    """
    This function finds an optimal placement for seismic scanners inside a pyramid,
    balancing the number of scans and the total volume covered.

    The strategy is a greedy algorithm executed in two phases:
    1.  Place large, long-range scanners (r=20m) first to cover as much volume as efficiently as possible.
    2.  Fill the remaining gaps with small, short-range scanners (r=7m).

    To make the search for valid, non-overlapping positions efficient, the algorithm uses:
    - A discrete grid for scanner coordinates (multiples of 0.5m), as specified.
    - A bottom-up, row-by-row search order (z, then y, then x) for deterministic placement.
    - A spatial hashing grid to speed up collision (overlap) checks, which is crucial for performance.
      A candidate sphere's position is checked only against existing spheres in its spatial vicinity,
      rather than the entire list of placed spheres.

    The final output is the number of long-range scans (n), short-range scans (m), and the
    percentage of the pyramid's volume covered by these scans (p), formatted as n:m:p.
    """
    
    # --- Parameters ---
    BASE_SIDE = 150.0  # Pyramid base side in meters
    HEIGHT = 110.0     # Pyramid height in meters
    R_LONG = 20.0      # Long-range scanner radius
    R_SHORT = 7.0      # Short-range scanner radius
    STEP = 0.5         # Coordinate grid step

    # --- Data structures for placed scans ---
    placed_scans = []
    
    # A spatial grid for fast overlap checking. Cell size is based on the larger scanner.
    # This choice works for both phases.
    SPATIAL_GRID_CELL_SIZE = 2 * R_LONG
    spatial_grid = {}

    def get_cell_index(point, cell_size):
        """Calculates the 3D index for a point in the spatial grid."""
        return (
            math.floor(point[0] / cell_size),
            math.floor(point[1] / cell_size),
            math.floor(point[2] / cell_size)
        )

    def is_contained(center, radius):
        """Checks if a sphere is fully contained within the pyramid."""
        cx, cy, cz = center
        
        # Check vertical bounds
        if not (cz - radius >= 0 and cz + radius <= HEIGHT):
            return False
        
        # Check horizontal bounds. The sphere must fit within the pyramid's cross-section
        # at the sphere's highest point (z = cz + radius), where the pyramid is narrowest.
        # This is the sufficient condition for containment within the side planes.
        pyramid_half_width_at_top = (BASE_SIDE / 2.0) * (1 - (cz + radius) / HEIGHT)
        if pyramid_half_width_at_top < 0: 
            return False
            
        if not (abs(cx) + radius <= pyramid_half_width_at_top and \
                abs(cy) + radius <= pyramid_half_width_at_top):
            return False
            
        return True

    def check_overlap_in_neighborhood(center, radius, existing_scans_grid, cell_size):
        """Checks for overlap with existing spheres in the local neighborhood."""
        center_cell_idx = get_cell_index(center, cell_size)
        
        # Iterate through the 3x3x3 block of cells around the candidate's cell
        for dz in range(-1, 2):
            for dy in range(-1, 2):
                for dx in range(-1, 2):
                    neighbor_cell_idx = (
                        center_cell_idx[0] + dx,
                        center_cell_idx[1] + dy,
                        center_cell_idx[2] + dz
                    )
                    
                    if neighbor_cell_idx in existing_scans_grid:
                        # Check against all scans in this neighboring cell
                        for neighbor_center, neighbor_radius in existing_scans_grid[neighbor_cell_idx]:
                            dist_sq = (center[0] - neighbor_center[0])**2 + \
                                      (center[1] - neighbor_center[1])**2 + \
                                      (center[2] - neighbor_center[2])**2
                            
                            radii_sum_sq = (radius + neighbor_radius)**2
                            if dist_sq < radii_sum_sq:
                                return True # Found an overlap
        return False # No overlap found

    def place_scans_of_radius(radius, current_scans, grid):
        """Greedy algorithm to place scans of a given radius."""
        
        # Define search bounds for the center coordinate z
        cz_min = radius
        cz_max = HEIGHT - radius
        
        # Iterate through possible z, y, x coordinates (bottom-up, deterministic order)
        z_coords = np.arange(cz_min, cz_max + STEP, STEP)
        for z in z_coords:
            # Determine max horizontal coordinate for a center at this height to narrow down the search
            max_abs_c = (BASE_SIDE / 2.0) * (1 - (z + radius) / HEIGHT) - radius
            if max_abs_c < 0:
                continue
            
            xy_coord_bound = math.floor(max_abs_c / STEP) * STEP
            xy_coords = np.arange(-xy_coord_bound, xy_coord_bound + STEP, STEP)

            for y in xy_coords:
                for x in xy_coords:
                    candidate_center = (x, y, z)
                    
                    # Using the loop bounds is a good first-pass check, but is_contained is fully robust.
                    if not is_contained(candidate_center, radius):
                        continue

                    if check_overlap_in_neighborhood(candidate_center, radius, grid, SPATIAL_GRID_CELL_SIZE):
                        continue

                    # If valid, place the scan and update data structures
                    scan_info = (candidate_center, radius)
                    current_scans.append(scan_info)
                    cell_idx = get_cell_index(candidate_center, SPATIAL_GRID_CELL_SIZE)
                    if cell_idx not in grid:
                        grid[cell_idx] = []
                    grid[cell_idx].append(scan_info)

    # --- Phase 1: Place Long-Range Scans ---
    place_scans_of_radius(R_LONG, placed_scans, spatial_grid)
    n_long = len(placed_scans)
    
    # --- Phase 2: Place Short-Range Scans ---
    place_scans_of_radius(R_SHORT, placed_scans, spatial_grid)
    n_short = len(placed_scans) - n_long

    # --- Calculate Final Results ---
    V_pyramid = (1.0/3.0) * BASE_SIDE**2 * HEIGHT
    V_long_single = (4.0/3.0) * math.pi * R_LONG**3
    V_short_single = (4.0/3.0) * math.pi * R_SHORT**3
    
    V_scanned_total = (n_long * V_long_single) + (n_short * V_short_single)
    
    coverage_ratio = (V_scanned_total / V_pyramid) * 100

    # --- Final Output ---
    # The problem requires printing the numbers for the final result.
    # Format: n:m:p
    # n: number of long range scans
    # m: number of short range scans
    # p: coverage ratio, rounded to 0.1
    final_n = n_long
    final_m = n_short
    final_p = round(coverage_ratio, 1)

    print(f"{final_n}:{final_m}:{final_p}")
    return f"{final_n}:{final_m}:{final_p}"

# --- Execution ---
if __name__ == "__main__":
    solve()