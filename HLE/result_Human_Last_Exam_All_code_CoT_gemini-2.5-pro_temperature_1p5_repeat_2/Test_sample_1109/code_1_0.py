import math
import numpy as np

# This script implements a greedy algorithm to solve the seismic scanning placement problem.
# It prioritizes placing larger, more efficient scanners first, then fills the
# remaining space with smaller scanners.
#
# The core logic involves:
# 1. Defining the pyramid's geometry and scanner properties.
# 2. A function to check if a spherical scan at a given location is valid, meaning:
#    a. It's entirely within the pyramid boundaries.
#    b. It doesn't overlap with any previously placed scans.
# 3. An optimized greedy placement loop that iterates through possible locations
#    (on a 0.5m grid) from the bottom of the pyramid upwards and from the center outwards.
#    This loop uses a spatial grid data structure to quickly check for collisions.
# 4. A final calculation of the total coverage volume and ratio.
#
# NOTE: This is a computationally intensive task. The execution may take a few minutes
# depending on the machine's performance.

def solve_seismic_scanning():
    """
    Finds an optimal placement of seismic scanners in the Isis pyramid.
    """
    
    # --- Pyramid and Scanner Parameters ---
    B_SIDE = 150.0  # Base side length
    HEIGHT = 110.0  # Height
    R_LONG = 20.0
    R_SHORT = 7.0
    GRID_STEP = 0.5

    # --- Pre-calculated constants for validity check ---
    # From the pyramid's slanted face plane equation: (2/B_SIDE)*x + (1/HEIGHT)*z = 1
    # which simplifies to 22x + 15z - 1650 = 0.
    # The distance from a point to this plane is used to check the boundary condition.
    SQRT_709 = math.sqrt(22**2 + 15**2)

    # --- Data structures ---
    placed_spheres = []
    # Use a spatial grid for fast collision detection. Cell size is 2 * R_LONG.
    cell_size = 2 * R_LONG
    grid = {}

    def get_cell_coords(x, y, z):
        """Calculates the grid cell indices for a given point."""
        # Offset coordinates to ensure non-negative indices for the dictionary key
        x_offset = B_SIDE
        y_offset = B_SIDE
        ix = int((x + x_offset) / cell_size)
        iy = int((y + y_offset) / cell_size)
        iz = int(z / cell_size) # z is always >= 0
        return (ix, iy, iz)

    def is_valid_placement(x, y, z, r):
        """Checks if a sphere can be placed at (x, y, z) with radius r."""
        # 1. Check if the sphere is geometrically inside the pyramid
        if not (r <= z <= HEIGHT - r):
            return False
        
        # This formula determines the max allowed |x| or |y| coordinate for a sphere's
        # center at height z, ensuring it stays within the slanted faces.
        max_abs_coord = (1650 - 15 * z - r * SQRT_709) / 22.0
        if abs(x) > max_abs_coord or abs(y) > max_abs_coord:
            return False
            
        # 2. Check for overlaps with already placed spheres using the spatial grid
        cell_coords = get_cell_coords(x, y, z)
        # Check the 3x3x3 cube of cells around the sphere's cell
        for di in range(-1, 2):
            for dj in range(-1, 2):
                for dk in range(-1, 2):
                    neighbor_cell = (cell_coords[0] + di, cell_coords[1] + dj, cell_coords[2] + dk)
                    if neighbor_cell in grid:
                        for p_sphere_idx in grid[neighbor_cell]:
                            px, py, pz, pr = placed_spheres[p_sphere_idx]
                            # Using squared distances to avoid costly sqrt operations
                            dist_sq = (x - px)**2 + (y - py)**2 + (z - pz)**2
                            min_dist_sq = (r + pr)**2
                            if dist_sq < min_dist_sq - 1e-9: # Epsilon for float precision
                                return False
        return True

    def place_sphere(x, y, z, r):
        """Adds a new sphere to our data structures."""
        sphere_idx = len(placed_spheres)
        placed_spheres.append((x, y, z, r))
        
        cell = get_cell_coords(x, y, z)
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(sphere_idx)

    def fill_pyramid_with_scans(r):
        """Iterates through grid points and places spheres of a given radius."""
        # Iterate from the lowest possible z to the highest
        z_min = math.ceil(r / GRID_STEP) * GRID_STEP
        z_max = math.floor((HEIGHT - r) / GRID_STEP) * GRID_STEP

        for z in np.arange(z_min, z_max + GRID_STEP, GRID_STEP):
            # Determine search bounds for this z-level to reduce iterations
            max_abs_coord = (1650 - 15 * z - r * SQRT_709) / 22.0
            if max_abs_coord < 0:
                continue

            coord_limit = math.floor(max_abs_coord / GRID_STEP) * GRID_STEP
            
            # Iterate from the center outwards to prioritize central placements
            y_coords = sorted(list(np.arange(-coord_limit, coord_limit + GRID_STEP, GRID_STEP)), key=abs)
            x_coords = sorted(list(np.arange(-coord_limit, coord_limit + GRID_STEP, GRID_STEP)), key=abs)
            
            for y in y_coords:
                for x in x_coords:
                    if is_valid_placement(x, y, z, r):
                        place_sphere(x, y, z, r)

    # --- Main Execution ---
    # Step 1: Place long-range scanners (r=20)
    fill_pyramid_with_scans(R_LONG)
    num_long_scans = len(placed_spheres)

    # Step 2: Fill remaining space with short-range scanners (r=7)
    fill_pyramid_with_scans(R_SHORT)
    num_short_scans = len(placed_spheres) - num_long_scans

    # Step 3: Calculate final results and coverage ratio
    V_pyramid = (1.0/3.0) * (B_SIDE**2) * HEIGHT
    V_long = (4.0/3.0) * math.pi * (R_LONG**3)
    V_short = (4.0/3.0) * math.pi * (R_SHORT**3)
    
    total_scanned_volume = (num_long_scans * V_long) + (num_short_scans * V_short)
    coverage_ratio = (total_scanned_volume / V_pyramid) * 100.0
    
    # Format and print the final output numbers
    n = num_long_scans
    m = num_short_scans
    p = round(coverage_ratio, 1)

    print("Optimal scanning configuration found:")
    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print(f"Total coverage ratio (p): {p}%")
    print("\nFinal Answer Format (n:m:p):")
    print(f"{n}:{m}:{p}")

# Run the simulation
solve_seismic_scanning()