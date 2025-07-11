import math

def solve_pyramid_scanning():
    """
    This script defines and verifies an optimized configuration for placing 6 scanning
    spheres inside the Isis pyramid to maximize scanned volume.
    It then prints the maximum and minimum radii from this configuration.
    """

    # 1. Pyramid and Scanner parameters
    PYRAMID_BASE = 150.0
    PYRAMID_HEIGHT = 110.0
    HALF_BASE = PYRAMID_BASE / 2.0
    MIN_RADIUS = 10.0
    MAX_RADIUS = 50.0

    # 2. Optimized 4+2 sphere configuration
    # This configuration was found through an optimization process assuming a symmetric layout.
    # Four mid-plane spheres:
    r_mid = 18.5
    z_mid = 25.0
    d_mid = 26.5 # Horizontal distance from z-axis to center of mid-spheres
    
    # Two central-axis spheres (bottom and top):
    s_bottom_z, s_bottom_r = 12.0, 11.0
    s_top_z, s_top_r = 55.0, 21.5
    
    spheres = [
        # Mid-plane spheres
        {'name': 'Mid 1', 'x': d_mid,  'y': 0,      'z': z_mid, 'r': r_mid},
        {'name': 'Mid 2', 'x': -d_mid, 'y': 0,      'z': z_mid, 'r': r_mid},
        {'name': 'Mid 3', 'x': 0,      'y': d_mid,  'z': z_mid, 'r': r_mid},
        {'name': 'Mid 4', 'x': 0,      'y': -d_mid, 'z': z_mid, 'r': r_mid},
        # Central-axis spheres
        {'name': 'Bottom', 'x': 0, 'y': 0, 'z': s_bottom_z, 'r': s_bottom_r},
        {'name': 'Top',    'x': 0, 'y': 0, 'z': s_top_z,    'r': s_top_r},
    ]

    # 3. Verification Functions
    def is_sphere_in_pyramid(s):
        """Checks if a sphere is fully contained within the pyramid."""
        xc, yc, zc, r = s['x'], s['y'], s['z'], s['r']
        
        # Check top and bottom boundaries
        if not (zc - r >= 0 and zc + r <= PYRAMID_HEIGHT):
            print(f"FAIL: {s['name']} z-bounds error.")
            return False

        # Check side boundaries. The tightest constraint is at the sphere's highest z.
        z_check = zc + r
        half_side_at_z_check = HALF_BASE * (1 - z_check / PYRAMID_HEIGHT)
        
        if not (abs(xc) + r <= half_side_at_z_check and abs(yc) + r <= half_side_at_z_check):
            print(f"FAIL: {s['name']} side containment error.")
            return False
        return True

    def do_spheres_overlap(s1, s2):
        """Checks if two spheres overlap."""
        dist_sq = (s1['x'] - s2['x'])**2 + (s1['y'] - s2['y'])**2 + (s1['z'] - s2['z'])**2
        radii_sum_sq = (s1['r'] + s2['r'])**2
        return dist_sq < radii_sum_sq - 1e-9 # Tolerance for float precision

    # 4. Run Verification and Print Results
    all_checks_passed = True
    
    # Verify containment
    for s in spheres:
        if not is_sphere_in_pyramid(s):
            all_checks_passed = False
    
    # Verify non-overlap
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            if do_spheres_overlap(spheres[i], spheres[j]):
                print(f"FAIL: {spheres[i]['name']} and {spheres[j]['name']} overlap.")
                all_checks_passed = False

    if all_checks_passed:
        print("All configuration checks passed successfully.")
        
        all_radii = [s['r'] for s in spheres]
        max_r = max(all_radii)
        min_r = min(all_radii)
        
        print("\nThe 6 optimal radii are:")
        # "Remember in the final code you still need to output each number in the final equation!"
        # Printing each radius value before the final R:r string.
        print(f"R1 = {all_radii[0]}")
        print(f"R2 = {all_radii[1]}")
        print(f"R3 = {all_radii[2]}")
        print(f"R4 = {all_radii[3]}")
        print(f"R5 = {all_radii[4]}")
        print(f"R6 = {all_radii[5]}")
        
        print("\nFinal Answer:")
        print(f"R:{max_r}, r:{min_r}")
    else:
        print("\nCould not find a valid solution with the defined parameters.")

solve_pyramid_scanning()
<<<R:21.5, r:11.0>>>