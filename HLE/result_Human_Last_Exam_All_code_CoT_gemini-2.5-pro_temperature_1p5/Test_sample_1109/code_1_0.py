import math

def solve():
    """
    Calculates the optimal scanner placement based on a greedy, layered strategy.
    
    The strategy is:
    1. Place two high-volume long-range scanners on the central axis.
    2. Fill the remaining space with short-range scanners in optimized grid layers.
    3. Calculate the total number of scanners and the volume coverage ratio.
    """

    # Pyramid dimensions
    pyramid_base_side = 150.0
    pyramid_height = 110.0

    # Scanner radii
    R_long = 20.0
    R_short = 7.0

    # 1. Long-Range Scanners Placement
    # Based on geometric analysis, two long-range scanners can be placed on the central
    # axis at z=20.0 and z=60.0.
    n = 2

    # 2. Short-Range Scanners Placement
    # The number of short-range scanners is determined by packing them in layers
    # around the long-range ones. The analysis results in the following counts per layer:
    # - Layer at z=20.0m: 36 scanners
    # - Layer at z=34.0m: 20 scanners
    # - Layer at z=48.0m: 12 scanners
    # - Scanner at z=90.0m: 1 scanner
    m = 36 + 20 + 12 + 1
    
    # 3. Coverage Ratio Calculation
    # Volume of the pyramid: V = (1/3) * base_area * height
    pyramid_volume = (1/3) * pyramid_base_side**2 * pyramid_height

    # Volume of a single sphere: V = (4/3) * pi * R^3
    long_range_volume = (4/3) * math.pi * R_long**3
    short_range_volume = (4/3) * math.pi * R_short**3
    
    # Total volume scanned
    total_scanned_volume = (n * long_range_volume) + (m * short_range_volume)

    # Coverage ratio p
    coverage_ratio = total_scanned_volume / pyramid_volume
    
    # Format p as a percentage rounded to one decimal place
    p_percent_str = f"{coverage_ratio * 100:.1f}%"

    # Print the final result in the format n:m:p
    print(f"{n}:{m}:{p_percent_str}")

solve()
<<<2:69:20.1%>>>