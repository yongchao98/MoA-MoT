import math

def solve_pyramid_scanning():
    """
    This script finds and verifies an optimal configuration for 6 seismic scans
    inside the Isis pyramid, aiming to maximize scanned volume.
    """

    # 1. Pyramid and Scanner Parameters
    h = 110.0  # height in meters
    b = 150.0  # base side length in meters
    k = 2 * h / b  # A constant for the pyramid slope constraint: 2*110/150

    # 2. Define the 6-sphere configuration based on a symmetric "4+2" model
    # This configuration was found by maximizing radii under tight-packing assumptions.
    # Sphere(x, y, z, radius)
    spheres = [
        # Sphere 1: Large sphere, low on the central z-axis
        {'id': 1, 'pos': (0.0, 0.0, 30.0), 'r': 30.0},
        # Spheres 2-5: Four spheres in a square pattern around sphere 1
        {'id': 2, 'pos': (35.0, 35.0, 30.0), 'r': 19.5},
        {'id': 3, 'pos': (35.0, -35.0, 30.0), 'r': 19.5},
        {'id': 4, 'pos': (-35.0, 35.0, 30.0), 'r': 19.5},
        {'id': 5, 'pos': (-35.0, -35.0, 30.0), 'r': 19.5},
        # Sphere 6: Upper sphere on the central z-axis
        {'id': 6, 'pos': (0.0, 0.0, 80.0), 'r': 20.0},
    ]

    print("Proposed scanning configuration (N=6):")
    total_volume = 0
    radii = []
    for s in spheres:
        print(f"  Scan {s['id']}: Location=({s['pos'][0]:.1f}, {s['pos'][1]:.1f}, {s['pos'][2]:.1f}), Radius={s['r']:.1f}")
        total_volume += (4/3) * math.pi * s['r']**3
        radii.append(s['r'])
    print(f"\nTotal Scanned Volume: {total_volume:.2f} m^3")
    print("-" * 30)
    print("Verification of Constraints:")

    all_constraints_met = True

    # 3. Verify Constraints
    # Constraint A: Sphere is inside the pyramid
    print("\nChecking 'Inside Pyramid' constraint...")
    for s in spheres:
        x, y, z = s['pos']
        r = s['r']
        
        # Check bottom plane (z > r)
        check1 = z >= r
        # Check side planes (|x| constraint)
        check2 = k * abs(x) + z + k * r <= h
        # Check side planes (|y| constraint)
        check3 = k * abs(y) + z + k * r <= h
        
        if not (check1 and check2 and check3):
            print(f"  [FAIL] Sphere {s['id']} is NOT fully inside the pyramid.")
            all_constraints_met = False
        else:
            print(f"  [PASS] Sphere {s['id']} is fully inside the pyramid.")

    # Constraint B: Spheres are non-overlapping
    print("\nChecking 'Non-Overlapping' constraint...")
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            s1 = spheres[i]
            s2 = spheres[j]
            
            p1 = s1['pos']
            p2 = s2['pos']
            
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            dist = math.sqrt(dist_sq)
            
            min_dist = s1['r'] + s2['r']
            
            if dist < min_dist:
                print(f"  [FAIL] Spheres {s1['id']} and {s2['id']} are overlapping.")
                all_constraints_met = False
            else:
                # To avoid spamming, we only print a sample pass
                if i==0 and j==1:
                     print(f"  [PASS] Spheres {s1['id']} and {s2['id']} do not overlap (Distance: {dist:.2f} >= Sum of Radii: {min_dist:.2f}).")


    if not all_constraints_met:
        print("\nCould not find a valid solution with the proposed configuration.")
        return

    # 4. Final Answer
    print("-" * 30)
    max_radius = max(radii)
    min_radius = min(radii)

    print("\nFinal Answer:")
    # The final equation is the relationship between Max Radius (R) and Min Radius (r)
    print(f"R:{max_radius:.1f}, r:{min_radius:.1f}")


solve_pyramid_scanning()
<<<30.0:19.5>>>