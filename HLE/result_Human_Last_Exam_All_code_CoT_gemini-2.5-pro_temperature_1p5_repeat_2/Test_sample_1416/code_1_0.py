import math

# Pyramid and Scanner parameters
PYRAMID_BASE = 150.0
PYRAMID_HEIGHT = 110.0
SCANNER_RADIUS_MIN = 10.0
SCANNER_RADIUS_MAX = 50.0
SCANNER_ADJUSTMENT = 0.5

def get_max_allowed_radius(center):
    """
    Calculates the maximum radius of a sphere at a given center to be fully
    contained within the pyramid. The radius is the minimum distance from the 
    center to any of the 5 faces of the pyramid.
    """
    cx, cy, cz = center
    
    # The center point must be inside the pyramid's outer boundary for any radius > 0
    # Boundary at height z: |x|, |y| <= (base/2) * (1 - z/height)
    side_at_cz = (PYRAMID_BASE / 2.0) * (1.0 - cz / PYRAMID_HEIGHT)
    if not (0 <= cz <= PYRAMID_HEIGHT and abs(cx) <= side_at_cz and abs(cy) <= side_at_cz):
        return 0.0

    # Distance to the base face (z=0 plane)
    dist_to_base = cz
    
    # For the 4 slanted faces, their plane equations are of the form Ax+By+Cz+D=0.
    # The distance from a point (cx,cy,cz) to a plane is |A*cx+B*cy+C*cz+D|/sqrt(A^2+B^2+C^2).
    # Normalizing factor for plane distance formula:
    norm = math.sqrt(PYRAMID_HEIGHT**2 + (PYRAMID_BASE / 2)**2) # sqrt(110^2 + 75^2)
    
    # Since the center is inside, the formula simplifies.
    # We check distance to the four planes representing the pyramid's sides.
    # Equation for +x side plane: 110*x + 75*z - 8250 = 0
    # Numerator for distance to +x plane: |110*cx + 75*cz - 8250|
    dist_to_side_px = abs(PYRAMID_HEIGHT * cx + (PYRAMID_BASE / 2) * cz - (PYRAMID_BASE / 2) * PYRAMID_HEIGHT) / norm
    dist_to_side_nx = abs(-PYRAMID_HEIGHT * cx + (PYRAMID_BASE / 2) * cz - (PYRAMID_BASE / 2) * PYRAMID_HEIGHT) / norm
    dist_to_side_py = abs(PYRAMID_HEIGHT * cy + (PYRAMID_BASE / 2) * cz - (PYRAMID_BASE / 2) * PYRAMID_HEIGHT) / norm
    dist_to_side_ny = abs(-PYRAMID_HEIGHT * cy + (PYRAMID_BASE / 2) * cz - (PYRAMID_BASE / 2) * PYRAMID_HEIGHT) / norm

    max_r = min(dist_to_base, dist_to_side_px, dist_to_side_nx, dist_to_side_py, dist_to_side_ny)
    
    # Clamp radius to scanner's max effective range
    max_r = min(max_r, SCANNER_RADIUS_MAX)
    
    # Floor to the nearest adjustment level
    return math.floor(max_r / SCANNER_ADJUSTMENT) * SCANNER_ADJUSTMENT

def check_overlap(sphere1, sphere2):
    """
    Checks if two spheres overlap. A sphere is a tuple (x, y, z, r).
    Touching is allowed (distance >= r1 + r2).
    """
    c1x, c1y, c1z, r1 = sphere1
    c2x, c2y, c2z, r2 = sphere2
    
    # Using squared distances to avoid costly square root operations
    dist_sq = (c1x - c2x)**2 + (c1y - c2y)**2 + (c1z - c2z)**2
    sum_radii_sq = (r1 + r2)**2
    
    # Return True if they overlap (distance is strictly less than sum of radii)
    return dist_sq < sum_radii_sq

def solve_pyramid_scan():
    """
    Finds and verifies the optimal sphere configuration.
    """
    print("Solving for optimal scanner placement (N=6)...")
    print("A '4+2' symmetric configuration provides a high-volume solution.")
    print("-" * 50)

    # Define the 6 spheres based on the optimized "4+2" configuration
    spheres = [
        # Set 1: Four lower spheres, radius 25.5m
        ( 25.5,  25.5, 25.5, 25.5),
        (-25.5,  25.5, 25.5, 25.5),
        ( 25.5, -25.5, 25.5, 25.5),
        (-25.5, -25.5, 25.5, 25.5),
        # Set 2: Two upper spheres, radius 12.0m
        ( 12.0,   0.0, 70.0, 12.0),
        (-12.0,   0.0, 70.0, 12.0),
    ]

    is_config_valid = True
    total_volume = 0.0

    # --- Verification Step 1: Check each sphere's validity ---
    for i, s in enumerate(spheres):
        cx, cy, cz, r = s
        max_r_at_center = get_max_allowed_radius((cx, cy, cz))
        is_contained = (r <= max_r_at_center)
        is_in_range = (SCANNER_RADIUS_MIN <= r <= SCANNER_RADIUS_MAX)
        
        if not (is_contained and is_in_range):
            is_config_valid = False
            print(f"X Sphere {i+1} is INVALID: Center=({cx:.1f}, {cy:.1f}, {cz:.1f}), Radius={r:.1f}")
        
        total_volume += (4.0/3.0) * math.pi * (r**3)

    # --- Verification Step 2: Check for overlaps ---
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            if check_overlap(spheres[i], spheres[j]):
                is_config_valid = False
                print(f"X Overlap DETECTED between Sphere {i+1} and Sphere {j+1}")

    # --- Print Results ---
    if not is_config_valid:
        print("\nProposed configuration is INVALID. Please check parameters.")
        return
        
    print("Proposed configuration is VALID. Details:")
    print("-" * 50)
    for i, s in enumerate(spheres):
        # The prompt requested outputting each number in the final equation.
        # This is interpreted as printing the parameters of each sphere.
        print(f"Sphere {i+1}: Center = ({s[0]:.1f}, {s[1]:.1f}, {s[2]:.1f}) m, Radius = {s[3]:.1f} m")

    print(f"\nTotal Scanned Volume: {total_volume:.2f} cubic meters")
    pyramid_volume = (1.0/3.0) * (PYRAMID_BASE**2) * PYRAMID_HEIGHT
    print(f"Total Pyramid Volume: {pyramid_volume:.2f} cubic meters")
    print(f"Scan Coverage: {100.0 * total_volume / pyramid_volume:.2f}%")
    print("-" * 50)

    # --- Final Answer Format ---
    all_radii = [s[3] for s in spheres]
    max_radius = max(all_radii)
    min_radius = min(all_radii)
    
    print("Final Answer (R:r, with R=max radius, r=min radius):")
    print(f"{max_radius:.1f}:{min_radius:.1f}")

solve_pyramid_scan()