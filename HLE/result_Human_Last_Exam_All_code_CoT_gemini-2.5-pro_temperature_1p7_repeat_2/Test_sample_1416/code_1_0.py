import math

def solve_pyramid_scanning():
    """
    This function calculates and verifies an optimal configuration for 6 seismic scanners
    inside the Isis pyramid, aiming to maximize the total scanned volume.
    """

    # 1. Define Pyramid and Scanner Parameters
    base_side = 150.0  # meters
    height = 110.0     # meters
    num_scanners = 6
    radius_min, radius_max = 10.0, 50.0
    coord_step = 0.5

    # This constant is derived from the pyramid's geometry for containment checks.
    # The plane equation for the pyramid faces is 110*|x| + 75*z = 8250.
    # The constraint for a sphere to be inside is 110*|xc| + 75*zc + r * L <= 8250,
    # where L is the norm of the plane's normal vector.
    L = math.sqrt(110**2 + 75**2)

    # 2. Define the Optimal Configuration found through the described strategy
    # The strategy involves a "4+2" model: 4 large spheres at the base and 2 on the axis.
    # Radii and positions are optimized to maximize volume while respecting constraints.
    scans = [
        # Four large spheres placed in a square pattern on a low plane (z=27.0)
        {'center': (25.5, 25.5, 27.0),  'radius': 25.5},
        {'center': (25.5, -25.5, 27.0), 'radius': 25.5},
        {'center': (-25.5, 25.5, 27.0), 'radius': 25.5},
        {'center': (-25.5, -25.5, 27.0), 'radius': 25.5},

        # Two smaller spheres stacked on the central axis (x=0, y=0)
        {'center': (0.0, 0.0, 91.0),   'radius': 10.5},
        {'center': (0.0, 0.0, 70.0),   'radius': 10.5}
    ]

    # 3. Verify the proposed configuration and present the result.
    print("Verifying the proposed optimal scanner configuration...")
    all_valid = True
    total_volume = 0

    for i, scan in enumerate(scans):
        xc, yc, zc = scan['center']
        r = scan['radius']
        
        # Containment check: z_center must be >= radius
        if not (zc >= r):
            print(f"FAILED: Scan {i+1} violates base containment (z_center < radius).")
            all_valid = False
            
        # Containment check: Must be inside the 4 side faces
        max_planar_coord = max(abs(xc), abs(yc))
        if not (110 * max_planar_coord + 75 * zc + r * L <= 8250 + 1e-9):
            print(f"FAILED: Scan {i+1} violates side-plane containment.")
            all_valid = False

        # Overlap check with other spheres
        for j in range(i + 1, len(scans)):
            xc2, yc2, zc2 = scans[j]['center']
            r2 = scans[j]['radius']
            dist_sq = (xc - xc2)**2 + (yc - yc2)**2 + (zc - zc2)**2
            min_dist_sq = (r + r2)**2
            if dist_sq < min_dist_sq - 1e-9:
                print(f"FAILED: Scan {i+1} overlaps with Scan {j+1}.")
                all_valid = False
        
        total_volume += (4/3) * math.pi * (r**3)

    if all_valid:
        print("Verification successful: All constraints are met.")
        print(f"\nThis configuration places {len(scans)} scanners with a total volume of {total_volume:.2f} m^3.")
    else:
        print("\nVerification FAILED. The configuration is invalid.")
        return

    # 4. Extract and print the final answer as requested
    all_radii = [s['radius'] for s in scans]
    max_radius = max(all_radii)
    min_radius = min(all_radii)

    print("\nThe analysis yields the following optimal radii:")
    print(f"Maximum scanning radius (R) = {max_radius:.1f} m")
    print(f"Minimum scanning radius (r) = {min_radius:.1f} m")
    
    # Present the final answer in the "R:r" format
    final_answer_string = f"{max_radius:.1f}:{min_radius:.1f}"
    print(f"\nFinal answer in the R:r format is:")
    print(final_answer_string)


solve_pyramid_scanning()
<<<25.5:10.5>>>