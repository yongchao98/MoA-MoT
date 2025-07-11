import math

# Define pyramid and scanner constants
PYRAMID_H = 110.0
PYRAMID_S = 150.0  # Base side length
R_LONG = 20.0
R_SHORT = 7.0
COORD_STEP = 0.5   # Scanner coordinates must be multiples of this value

def get_pyramid_volume():
    """Calculates the volume of the pyramid."""
    return (1/3) * PYRAMID_S**2 * PYRAMID_H

def get_sphere_volume(radius):
    """Calculates the volume of a sphere."""
    return (4/3) * math.pi * radius**3

def is_valid_in_pyramid(center, radius):
    """
    Checks if a sphere is completely inside the pyramid using the geometric constraints.
    A sphere is inside if its center is at least `radius` distance from all 5 faces.
    """
    x, y, z = center
    
    # 1. Check distance from base (z=0) and top (z=H)
    if not (radius <= z <= PYRAMID_H - radius):
        return False
        
    # 2. Check distance from the 4 slanted faces.
    # The condition is derived from the plane equation of each face.
    # For the +x face, the plane is H*x + (S/2)*z - H*S/2 = 0.
    # The distance from the center to this plane must be >= radius.
    
    # Length of the normal vector to the slanted face plane
    norm_len = math.sqrt(PYRAMID_H**2 + (PYRAMID_S/2.0)**2)
    # The condition simplifies to: H*|coord| + (S/2)*z <= H*S/2 - r*norm_len
    rhs = PYRAMID_H * (PYRAMID_S / 2.0) - radius * norm_len

    if PYRAMID_H * abs(x) + (PYRAMID_S / 2.0) * z > rhs:
        return False
    if PYRAMID_H * abs(y) + (PYRAMID_S / 2.0) * z > rhs:
        return False
        
    return True

def find_optimal_scans():
    """
    Finds an optimal set of scan locations using a greedy placement algorithm.
    """
    placed_scanners = []
    
    # Define scan types to try, large radius first for volume efficiency.
    # A search_step is used to create a grid of candidate points to check.
    scan_configs = [
        {"radius": R_LONG, "search_step": 2.0},
        {"radius": R_SHORT, "search_step": 1.0},
    ]

    for config in scan_configs:
        radius = config["radius"]
        search_step = config["search_step"]
        
        # Iterate through a grid of potential center points (z, y, x).
        # We search only in the positive quadrant (x>=0, y>=0) and use symmetry.
        z = radius
        while z <= PYRAMID_H - radius:
            # Max coordinate for the center at this height z
            # This is an optimization to avoid searching outside the pyramid.
            # A simplified boundary is used here. The strict check is in is_valid_in_pyramid.
            max_coord_at_z = (PYRAMID_S / 2.0) * (PYRAMID_H - z) / PYRAMID_H - radius
            if max_coord_at_z < 0:
                z += search_step
                continue

            y = 0.0
            while y <= max_coord_at_z:
                x = 0.0
                while x <= max_coord_at_z:
                    # Generate symmetric points (±x, ±y, z)
                    signs_x = [1] if x == 0 else [1, -1]
                    signs_y = [1] if y == 0 else [1, -1]
                    
                    for sx in signs_x:
                        for sy in signs_y:
                            px, py, pz = x * sx, y * sy, z
                            
                            # Round center to the nearest multiple of COORD_STEP (0.5m)
                            center = (
                                round(px / COORD_STEP) * COORD_STEP,
                                round(py / COORD_STEP) * COORD_STEP,
                                round(pz / COORD_STEP) * COORD_STEP
                            )
                            
                            if not is_valid_in_pyramid(center, radius):
                                continue

                            # Check for overlap with already placed scanners.
                            # Two spheres are non-overlapped if distance between centers >= sum of radii.
                            is_overlapping = False
                            for p_center, p_radius in placed_scanners:
                                dist_sq = (center[0] - p_center[0])**2 + (center[1] - p_center[1])**2 + (center[2] - p_center[2])**2
                                min_dist_sq = (radius + p_radius)**2
                                # Use a small tolerance for floating point comparisons
                                if dist_sq < min_dist_sq - 1e-9:
                                    is_overlapping = True
                                    break
                            
                            if not is_overlapping:
                                placed_scanners.append((center, radius))
                    x += search_step
                y += search_step
            z += search_step
            
    # Count the number of scans for each mode
    n_long = sum(1 for _, r in placed_scanners if r == R_LONG)
    n_short = sum(1 for _, r in placed_scanners if r == R_SHORT)

    # Calculate volumes and coverage ratio
    vol_long_single = get_sphere_volume(R_LONG)
    vol_short_single = get_sphere_volume(R_SHORT)
    total_scanned_volume = n_long * vol_long_single + n_short * vol_short_single
    pyramid_volume = get_pyramid_volume()
    coverage_ratio = (total_scanned_volume / pyramid_volume) * 100
    
    # Print the detailed breakdown of the coverage calculation
    print("Coverage Calculation:")
    print(f"Pyramid Volume = {pyramid_volume:.2f} m^3")
    print(f"Long-Range Scan Volume (each) = {vol_long_single:.2f} m^3")
    print(f"Short-Range Scan Volume (each) = {vol_short_single:.2f} m^3")
    print(f"Total Scanned Volume = ( {n_long} * {vol_long_single:.2f} ) + ( {n_short} * {vol_short_single:.2f} ) = {total_scanned_volume:.2f} m^3")
    print(f"Coverage Ratio = ({total_scanned_volume:.2f} / {pyramid_volume:.2f}) * 100% = {coverage_ratio:.1f}%")
    print("-" * 20)
    
    # Print the final result in the requested format
    print("Final Answer:")
    print(f"{n_long}:{n_short}:{coverage_ratio:.1f}")
    
    # Return the value for the final answer block
    return f"{n_long}:{n_short}:{coverage_ratio:.1f}"

if __name__ == '__main__':
    final_answer = find_optimal_scans()
    print(f"\n<<< {final_answer} >>>")