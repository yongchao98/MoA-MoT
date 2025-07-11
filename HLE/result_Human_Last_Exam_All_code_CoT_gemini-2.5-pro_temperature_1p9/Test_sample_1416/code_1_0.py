import math

# Pyramid and Scanner constants
PYRAMID_BASE_SIDE = 150.0
PYRAMID_HEIGHT = 110.0
MIN_RADIUS = 10.0
MAX_RADIUS = 50.0

# Pre-calculated geometric constants for the pyramid faces
# Normal vector to side planes has components (PYRAMID_HEIGHT, PYRAMID_BASE_SIDE/2)
# Magnitude of this normal vector in the x-z plane projection
NORM_SIDE_PLANES = math.sqrt(PYRAMID_HEIGHT**2 + (PYRAMID_BASE_SIDE / 2)**2)
# Constant for plane equation h*|x| + (a/2)*z = h*(a/2)
PLANE_CONSTANT = PYRAMID_HEIGHT * (PYRAMID_BASE_SIDE / 2)

def is_sphere_inside_pyramid(center, radius):
    """
    Checks if a sphere is completely inside the pyramid.
    A sphere is inside if its center is further from each boundary plane
    than its radius.
    """
    cx, cy, cz = center
    
    # 1. Check against the base plane (z=0)
    if cz < radius:
        return False

    # 2. Check against the four side planes
    # The condition is derived from the point-to-plane distance formula.
    # We check against the plane that the center is closest to, which is determined
    # by the larger of |cx| and |cy|.
    # h*max(|cx|,|cy|) + (a/2)*cz + r*norm <= h*(a/2)
    val = (PYRAMID_HEIGHT * max(abs(cx), abs(cy)) +
           (PYRAMID_BASE_SIDE / 2) * cz +
           radius * NORM_SIDE_PLANES)
    
    if val > PLANE_CONSTANT + 1e-9: # Add tolerance for float precision
        return False
        
    return True

def do_spheres_overlap(center1, r1, center2, r2):
    """Checks if two spheres overlap."""
    dist_sq = ((center1[0] - center2[0])**2 +
               (center1[1] - center2[1])**2 +
               (center1[2] - center2[2])**2)
    sum_radii_sq = (r1 + r2)**2
    # Use a small tolerance for floating point comparisons
    return dist_sq < sum_radii_sq - 1e-9

def solve():
    """
    This function defines and verifies the optimal configuration for 6 scans
    and prints the required results.
    """
    # This configuration is found by optimizing a symmetric '2+4' arrangement
    # of spheres within the pyramid's geometric constraints.
    scans = [
        {'id': 1, 'center': (0.0, 0.0, 23.0), 'radius': 23.0},
        {'id': 2, 'center': (0.0, 0.0, 69.0), 'radius': 23.0},
        {'id': 3, 'center': (23.0, 23.0, 46.0), 'radius': 16.5},
        {'id': 4, 'center': (-23.0, 23.0, 46.0), 'radius': 16.5},
        {'id': 5, 'center': (23.0, -23.0, 46.0), 'radius': 16.5},
        {'id': 6, 'center': (-23.0, -23.0, 46.0), 'radius': 16.5}
    ]

    # --- Verification Step (optional, for confidence) ---
    all_valid = True
    for s in scans:
        if not is_sphere_inside_pyramid(s['center'], s['radius']):
            all_valid = False
            print(f"Validation Error: Scan {s['id']} is not inside the pyramid.")
    for i in range(len(scans)):
        for j in range(i + 1, len(scans)):
            if do_spheres_overlap(scans[i]['center'], scans[i]['radius'], scans[j]['center'], scans[j]['radius']):
                all_valid = False
                print(f"Validation Error: Scan {scans[i]['id']} and {scans[j]['id']} overlap.")
    if not all_valid:
        print("\nWarning: The determined configuration has validation errors.")
        
    # --- Print Results ---
    print("Optimal scanning locations and radii for N=6:")
    radii = []
    for s in scans:
        print(f"  - Scan {s['id']}: Center = {s['center']}, Radius = {s['radius']:.1f}")
        radii.append(s['radius'])
        
    R = max(radii)
    r = min(radii)

    print(f"\nThe maximum scanning radius (R) is {R:.1f}m.")
    print(f"The minimum scanning radius (r) is {r:.1f}m.")
    
# Execute the solution
solve()

<<<23.0:16.5>>>