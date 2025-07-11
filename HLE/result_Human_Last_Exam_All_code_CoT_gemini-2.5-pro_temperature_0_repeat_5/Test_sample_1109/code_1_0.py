import math

def solve_seismic_scanning():
    """
    Calculates an optimal scanning configuration for the Isis pyramid.
    """
    # 1. Define constants for the pyramid and scanners
    PYRAMID_HEIGHT = 110.0
    PYRAMID_BASE_SIDE = 150.0
    PYRAMID_HALF_BASE = PYRAMID_BASE_SIDE / 2.0
    
    LONG_RANGE_RADIUS = 20.0
    SHORT_RANGE_RADIUS = 7.0

    # 2. Calculate total volume of the pyramid and individual scan volumes
    PYRAMID_VOLUME = (PYRAMID_BASE_SIDE**2 * PYRAMID_HEIGHT) / 3.0
    LONG_RANGE_VOLUME = (4.0/3.0) * math.pi * (LONG_RANGE_RADIUS**3)
    SHORT_RANGE_VOLUME = (4.0/3.0) * math.pi * (SHORT_RANGE_RADIUS**3)

    # 3. Determine the number and placement of scanners based on a structured packing strategy

    # --- Long-Range Scans (n) ---
    # Place two spheres on the central axis, pushed as high as possible to leave
    # a large contiguous volume at the wider base of the pyramid.
    # The centers are at (0, 0, 34.5) and (0, 0, 74.5).
    n = 2
    center_large_1_z = 34.5
    
    # --- Short-Range Scans (m) ---
    # Place small spheres in horizontal layers, on a square grid.
    m = 0
    
    # Geometric factor to determine the 'safe' inner pyramid for sphere centers
    PYRAMID_GEO_FACTOR = math.sqrt(1 + (PYRAMID_HALF_BASE / PYRAMID_HEIGHT)**2)
    GRID_SPACING = 2 * SHORT_RANGE_RADIUS  # 14.0m

    # Layer 1: z = 7.5m
    z1 = 7.5
    # Max horizontal coordinate for a center at this height
    max_x1 = PYRAMID_HALF_BASE * (1 - z1 / PYRAMID_HEIGHT) - SHORT_RANGE_RADIUS * PYRAMID_GEO_FACTOR
    # Max grid index (e.g., i in (i * 14, j * 14, z))
    max_i1 = int(max_x1 / GRID_SPACING)
    # Number of spheres in a (2*i+1) by (2*i+1) grid
    m_layer1 = (2 * max_i1 + 1)**2
    m += m_layer1

    # Layer 2: z = 21.5m (7.5 + 14.0)
    z2 = 21.5
    max_x2 = PYRAMID_HALF_BASE * (1 - z2 / PYRAMID_HEIGHT) - SHORT_RANGE_RADIUS * PYRAMID_GEO_FACTOR
    max_i2 = int(max_x2 / GRID_SPACING)
    # This layer has a "hole" due to the first large sphere. We subtract the 9 central grid points.
    m_layer2 = (2 * max_i2 + 1)**2 - 9
    m += m_layer2

    # Layer 3: z = 35.5m (21.5 + 14.0)
    z3 = 35.5
    max_x3 = PYRAMID_HALF_BASE * (1 - z3 / PYRAMID_HEIGHT) - SHORT_RANGE_RADIUS * PYRAMID_GEO_FACTOR
    max_i3 = int(max_x3 / GRID_SPACING)
    # This layer also has a hole. We subtract the 9 central grid points.
    m_layer3 = (2 * max_i3 + 1)**2 - 9
    m += m_layer3

    # Layer 4: z = 49.5m (35.5 + 14.0)
    z4 = 49.5
    max_x4 = PYRAMID_HALF_BASE * (1 - z4 / PYRAMID_HEIGHT) - SHORT_RANGE_RADIUS * PYRAMID_GEO_FACTOR
    max_i4 = int(max_x4 / GRID_SPACING)
    # This layer also has a hole. We subtract the 9 central grid points.
    m_layer4 = (2 * max_i4 + 1)**2 - 9
    m += m_layer4

    # Higher layers do not fit any spheres using this grid-based strategy.

    # 4. Calculate total coverage and format the output
    total_scanned_volume = (n * LONG_RANGE_VOLUME) + (m * SHORT_RANGE_VOLUME)
    coverage_ratio = (total_scanned_volume / PYRAMID_VOLUME) * 100
    
    # Round coverage to one decimal place as requested ("rounded to 0.1%")
    p_rounded = round(coverage_ratio, 1)

    # Print the final result in the format n:m:p
    print(f"{n}:{m}:{p_rounded}")

solve_seismic_scanning()
<<<2:177:38.9>>>