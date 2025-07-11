import math

def solve_seismic_scanning():
    """
    Calculates an optimal set of seismic scans for the Isis pyramid
    using a greedy placement algorithm.
    """
    # --- Pyramid and Scanner Parameters ---
    H = 110.0
    B_HALF = 75.0
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5
    
    # Use a spatial grid for efficient overlap detection.
    # Cell size is based on the largest radius to ensure effectiveness.
    CELL_SIZE = 2 * R_LONG 
    
    placed_spheres = []
    spatial_grid = {}

    def get_pyramid_half_side(z):
        """Calculates the half-side length of the pyramid's square cross-section at height z."""
        if not (0 <= z <= H):
            return 0
        return B_HALF * (1 - z / H)

    def is_sphere_inside(cx, cy, cz, r):
        """
        Conservatively checks if a sphere is entirely inside the pyramid
        by ensuring its axis-aligned bounding box fits.
        """
        # Check vertical bounds
        if not (r <= cz <= H - r):
            return False
        
        # The most restrictive horizontal check is at the sphere's lowest point.
        pyramid_side_at_bottom = get_pyramid_half_side(cz - r)
        if abs(cx) + r > pyramid_side_at_bottom:
            return False
        if abs(cy) + r > pyramid_side_at_bottom:
            return False
        return True

    def get_cell_indices(x, y, z):
        """Calculates the spatial grid cell indices for a given coordinate."""
        return (
            math.floor(x / CELL_SIZE),
            math.floor(y / CELL_SIZE),
            math.floor(z / CELL_SIZE)
        )

    def add_sphere_to_grid(center, radius):
        """Adds a newly placed sphere to the main list and the spatial grid."""
        sphere_data = (center, radius)
        placed_spheres.append(sphere_data)
        cell_indices = get_cell_indices(*center)
        if cell_indices not in spatial_grid:
            spatial_grid[cell_indices] = []
        spatial_grid[cell_indices].append(sphere_data)
        
    def check_for_overlaps(cx, cy, cz, r):
        """
        Checks if a candidate sphere overlaps with any previously placed spheres
        by querying the spatial grid for nearby objects.
        """
        center_cell = get_cell_indices(cx, cy, cz)
        
        # Check spheres in the 3x3x3 neighborhood of cells around the candidate.
        for ix_offset in range(-1, 2):
            for iy_offset in range(-1, 2):
                for iz_offset in range(-1, 2):
                    cell_to_check = (
                        center_cell[0] + ix_offset,
                        center_cell[1] + iy_offset,
                        center_cell[2] + iz_offset,
                    )
                    
                    if cell_to_check in spatial_grid:
                        for p_center, p_radius in spatial_grid[cell_to_check]:
                            dist_sq = (
                                (cx - p_center[0])**2 +
                                (cy - p_center[1])**2 +
                                (cz - p_center[2])**2
                            )
                            min_dist_sq = (r + p_radius)**2
                            if dist_sq < min_dist_sq:
                                return True # Overlap found
        return False # No overlap

    def place_scans(radius):
        """
        The main greedy placement function. It iterates through all valid grid 
        points and places a sphere if it fits and doesn't overlap.
        """
        # Define the search range for sphere centers on the 0.5m grid.
        z_min_grid = math.ceil(radius / GRID_STEP)
        z_max_grid = math.floor((H - radius) / GRID_STEP)
        
        # Iterate z from bottom-up
        for z_grid in range(z_min_grid, z_max_grid + 1):
            cz = z_grid * GRID_STEP
            
            # Determine max coordinate for a center at this height.
            max_abs_coord = get_pyramid_half_side(cz - radius) - radius
            if max_abs_coord < 0:
                continue
            
            xy_max_grid = math.floor(max_abs_coord / GRID_STEP)
            
            # Iterate x and y in a raster scan pattern.
            for x_grid in range(-xy_max_grid, xy_max_grid + 1):
                cx = x_grid * GRID_STEP
                for y_grid in range(-xy_max_grid, xy_max_grid + 1):
                    cy = y_grid * GRID_STEP
                    
                    # Optimization: rough check before detailed ones.
                    if abs(cx) > max_abs_coord or abs(cy) > max_abs_coord:
                        continue
                        
                    # 1. Perform detailed check that sphere is inside the pyramid.
                    if not is_sphere_inside(cx, cy, cz, radius):
                        continue

                    # 2. Check for overlaps with already placed spheres.
                    if check_for_overlaps(cx, cy, cz, radius):
                        continue
                        
                    # If all checks pass, place the sphere.
                    add_sphere_to_grid((cx, cy, cz), radius)

    # --- Execution ---
    # Phase 1: Place long-range scanners first for maximum volume efficiency.
    place_scans(R_LONG)
    num_long_range = len(placed_spheres)

    # Phase 2: Place short-range scanners to fill in the gaps.
    place_scans(R_SHORT)
    num_short_range = len(placed_spheres) - num_long_range
    
    # --- Final Calculation and Output ---
    volume_pyramid = (1/3) * (B_HALF * 2)**2 * H
    
    volume_scanned = 0
    for center, radius in placed_spheres:
        volume_scanned += (4/3) * math.pi * radius**3
        
    coverage_ratio = (volume_scanned / volume_pyramid) * 100
    
    n = num_long_range
    m = num_short_range
    p = round(coverage_ratio, 1)

    print(f"{n}:{m}:{p}")

# Run the simulation. Note: This is a computationally intensive task.
solve_seismic_scanning()
<<<10:88:55.3>>>