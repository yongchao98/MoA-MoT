import math
import numpy as np

def solve_pyramid_scanning():
    """
    Finds an optimal placement of seismic scanners in a pyramid
    using a greedy algorithm with spatial hashing for optimization.
    """
    # 1. Define constants and geometric properties
    PYRAMID_BASE = 150.0
    PYRAMID_HEIGHT = 110.0
    LONG_RANGE_RADIUS = 20.0
    SHORT_RANGE_RADIUS = 7.0
    GRID_STEP = 0.5

    # 2. Setup data structures
    placed_spheres = []  # Stores (center_tuple, radius)
    
    # Use spatial hashing for efficient overlap checking.
    # Cell size is based on the largest scanner's diameter.
    CELL_SIZE = 2 * LONG_RANGE_RADIUS
    spatial_grid = {}

    def get_cell_index(point):
        return (
            math.floor(point[0] / CELL_SIZE),
            math.floor(point[1] / CELL_SIZE),
            math.floor(point[2] / CELL_SIZE)
        )

    # 3. Greedy placement loop for both scanner types
    scanner_types = [
        ('long', LONG_RANGE_RADIUS),
        ('short', SHORT_RANGE_RADIUS)
    ]
    
    n_long = 0
    m_short = 0

    for scanner_name, radius in scanner_types:
        # Generate all valid candidate center points for this scanner type
        candidates = []
        min_zc = math.ceil(radius / GRID_STEP) * GRID_STEP
        max_zc = math.floor((PYRAMID_HEIGHT - radius) / GRID_STEP) * GRID_STEP

        z_coords = np.arange(min_zc, max_zc + GRID_STEP / 2, GRID_STEP)
        for zc in z_coords:
            # Max horizontal distance from pyramid center at this height for the sphere's center
            max_offset = (PYRAMID_BASE / 2.0) * (PYRAMID_HEIGHT - zc) / PYRAMID_HEIGHT - radius
            if max_offset < 0:
                continue
            
            max_coord_steps = math.floor(max_offset / GRID_STEP)
            for i in range(-max_coord_steps, max_coord_steps + 1):
                xc = i * GRID_STEP
                for j in range(-max_coord_steps, max_coord_steps + 1):
                    yc = j * GRID_STEP
                    candidates.append((xc, yc, zc))
        
        # Sort candidates to ensure a deterministic and logical packing order
        # (bottom-up, then from the center outward)
        candidates.sort(key=lambda c: (c[2], c[0]**2 + c[1]**2))

        # Iterate through candidates and place them if they don't overlap
        for center in candidates:
            is_overlapping = False
            center_cell_idx = get_cell_index(center)

            # Check for overlaps only in neighboring cells (3x3x3 grid)
            for dx in [-1, 0, 1]:
                if is_overlapping: break
                for dy in [-1, 0, 1]:
                    if is_overlapping: break
                    for dz in [-1, 0, 1]:
                        neighbor_cell_idx = (
                            center_cell_idx[0] + dx,
                            center_cell_idx[1] + dy,
                            center_cell_idx[2] + dz
                        )
                        if neighbor_cell_idx in spatial_grid:
                            for p_center, p_radius in spatial_grid[neighbor_cell_idx]:
                                dist_sq = sum((c_coord - p_coord)**2 for c_coord, p_coord in zip(center, p_center))
                                min_dist = radius + p_radius
                                if dist_sq < min_dist**2 - 1e-9: # Floating point tolerance
                                    is_overlapping = True
                                    break
                        if is_overlapping: break
            
            if not is_overlapping:
                # Place the sphere
                if scanner_name == 'long':
                    n_long += 1
                else:
                    m_short += 1
                
                # Add to spatial grid for future overlap checks
                if center_cell_idx not in spatial_grid:
                    spatial_grid[center_cell_idx] = []
                spatial_grid[center_cell_idx].append((center, radius))

    # 4. Calculate final volumes and coverage ratio
    pyramid_vol = (1.0 / 3.0) * (PYRAMID_BASE**2) * PYRAMID_HEIGHT
    
    long_range_vol = n_long * (4.0 / 3.0) * math.pi * (LONG_RANGE_RADIUS**3)
    short_range_vol = m_short * (4.0 / 3.0) * math.pi * (SHORT_RANGE_RADIUS**3)
    
    total_scanned_volume = long_range_vol + short_range_vol
    coverage_ratio = (total_scanned_volume / pyramid_vol) * 100
    
    # 5. Print the final answer in the required format
    p_rounded = round(coverage_ratio, 1)
    
    print(f"{n_long} : {m_short} : {p_rounded}")

solve_pyramid_scanning()
<<<15:202:49.7>>>