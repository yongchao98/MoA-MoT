import numpy as np
import math

def solve_seismic_scanning():
    """
    Solves the seismic scanning puzzle by simulating the placement of scanners
    inside a pyramid.
    """
    # 1. Define constants and pyramid/scanner properties
    H = 110.0
    BASE_SIDE = 150.0
    BASE_HALF_SIDE = BASE_SIDE / 2.0
    R_L = 20.0
    R_S = 7.0
    GRID_STEP = 0.5

    # --- Volume Calculations ---
    pyramid_vol = (1.0 / 3.0) * (BASE_SIDE**2) * H
    sphere_vol_L = (4.0 / 3.0) * math.pi * (R_L**3)
    sphere_vol_S = (4.0 / 3.0) * math.pi * (R_S**3)

    # 2. Place Long-Range Scanners (n)
    # This configuration is pre-determined to be an efficient packing.
    placed_long_range = []
    long_range_centers = []
    # Layer 1: z=20m, 3x3 grid with 40m spacing
    for i in [-1.0, 0.0, 1.0]:
        for j in [-1.0, 0.0, 1.0]:
            long_range_centers.append((i * 40.0, j * 40.0, 20.0))
    # Layer 2: z=60m, central scanner
    long_range_centers.append((0.0, 0.0, 60.0))
    
    for center in long_range_centers:
        placed_long_range.append({'center': center, 'radius': R_L})
    
    n = len(placed_long_range)
    scanned_vol = n * sphere_vol_L

    # 3. Prepare for Short-Range Scanner Placement
    # Create a 3D boolean grid representing all possible center locations for short-range scanners.
    
    # Grid boundaries
    z_min_grid, z_max_grid = R_S, H - R_S
    # The grid covers the full base to simplify indexing
    x_min_grid, x_max_grid = -BASE_HALF_SIDE, BASE_HALF_SIDE
    y_min_grid, y_max_grid = -BASE_HALF_SIDE, BASE_HALF_SIDE

    # Grid dimensions
    nz = int(round((z_max_grid - z_min_grid) / GRID_STEP)) + 1
    nx = int(round((x_max_grid - x_min_grid) / GRID_STEP)) + 1
    ny = int(round((y_max_grid - y_min_grid) / GRID_STEP)) + 1
    
    is_available = np.full((nx, ny, nz), True, dtype=bool)

    # Helper functions to map between world coordinates and grid indices
    def coord_to_idx(p):
        ix = int(round((p[0] - x_min_grid) / GRID_STEP))
        iy = int(round((p[1] - y_min_grid) / GRID_STEP))
        iz = int(round((p[2] - z_min_grid) / GRID_STEP))
        return (ix, iy, iz)

    def idx_to_coord(idx):
        x = x_min_grid + idx[0] * GRID_STEP
        y = y_min_grid + idx[1] * GRID_STEP
        z = z_min_grid + idx[2] * GRID_STEP
        return (x, y, z)

    # --- Block invalid grid points ---
    
    # Block points outside the pyramid or too close to long-range scanners
    print("Initializing grid and blocking invalid zones...")
    for iz in range(nz):
        z = z_min_grid + iz * GRID_STEP
        half_side_at_z = BASE_HALF_SIDE * (H - z) / H
        max_coord_at_z = half_side_at_z - R_S
        
        for ix in range(nx):
            x = x_min_grid + ix * GRID_STEP
            if abs(x) > max_coord_at_z:
                is_available[ix, :, iz] = False
                continue
            
            for iy in range(ny):
                if not is_available[ix, iy, iz]:
                    continue
                y = y_min_grid + iy * GRID_STEP
                if abs(y) > max_coord_at_z:
                    is_available[ix, iy, iz] = False
                    continue
                
                # Check overlap with long-range scanners
                for lr_sphere in placed_long_range:
                    dist_sq = (x - lr_sphere['center'][0])**2 + \
                              (y - lr_sphere['center'][1])**2 + \
                              (z - lr_sphere['center'][2])**2
                    if dist_sq < (R_S + R_L)**2:
                        is_available[ix, iy, iz] = False
                        break

    # 4. Greedy Placement of Short-Range Scanners (m)
    print("Placing short-range scanners...")
    m = 0
    block_radius_idx = int(math.ceil((2 * R_S) / GRID_STEP))

    # Iterate through the grid from bottom up
    for iz in range(nz):
        for ix in range(nx):
            for iy in range(ny):
                if is_available[ix, iy, iz]:
                    # Found a valid spot, place a sphere
                    m += 1
                    
                    # Block out nearby locations
                    center_coord = idx_to_coord((ix, iy, iz))
                    
                    min_ix = max(0, ix - block_radius_idx)
                    max_ix = min(nx, ix + block_radius_idx + 1)
                    min_iy = max(0, iy - block_radius_idx)
                    max_iy = min(ny, iy + block_radius_idx + 1)
                    min_iz = max(0, iz - block_radius_idx)
                    max_iz = min(nz, iz + block_radius_idx + 1)

                    for jz in range(min_iz, max_iz):
                        for jx in range(min_ix, max_ix):
                            for jy in range(min_iy, max_iy):
                                if is_available[jx, jy, jz]:
                                    p_coord = idx_to_coord((jx, jy, jz))
                                    dist_sq = (p_coord[0] - center_coord[0])**2 + \
                                              (p_coord[1] - center_coord[1])**2 + \
                                              (p_coord[2] - center_coord[2])**2
                                    if dist_sq < (2 * R_S)**2:
                                        is_available[jx, jy, jz] = False
    
    scanned_vol += m * sphere_vol_S
    coverage_ratio = (scanned_vol / pyramid_vol) * 100
    
    # 5. Output the results
    print("\n--- Scan Optimization Results ---")
    print(f"Number of long-range scans (n): {n}")
    print(f"Number of short-range scans (m): {m}")
    print(f"Total scanned volume: {scanned_vol:.1f} m^3")
    print(f"Total pyramid volume: {pyramid_vol:.1f} m^3")
    print(f"Coverage ratio (p): {coverage_ratio:.1f}%")
    
    # Final answer in the required format
    print("\nFinal Answer:")
    print(f"{n}:{m}:{coverage_ratio:.1f}")

solve_seismic_scanning()