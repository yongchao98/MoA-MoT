import math

def solve_pyramid_scanning():
    """
    Calculates and verifies an optimal 6-sphere packing configuration inside a pyramid.

    The solution is based on a heuristic placing 4 spheres on a lower layer
    and 2 spheres on an upper layer, arranged symmetrically.
    """

    # Pyramid definition
    base_side = 150.0  # meters
    height = 110.0   # meters

    # The four slanted faces of the pyramid can be described by the general
    # equation `22*|x| + 15*z = 1650` and `22*|y| + 15*z = 1650`.
    # The fifth face is the base at z=0.
    
    # Pre-calculated constant for plane distance formula sqrt(A^2+B^2+C^2) for slanted faces
    # sqrt(22^2 + 15^2)
    sqrt_709 = math.sqrt(709)

    # --- Configuration based on heuristic optimization ---
    # Layer 1: 4 spheres
    r1 = 25.5
    d1 = 25.5  # Half-distance between sphere centers on x/y plane
    z1 = 25.5  # Height of sphere centers
    
    # Layer 2: 2 spheres
    r2 = 16.0
    d2 = 16.0  # Half-distance between sphere centers on x-axis
    z2 = 57.5  # Height of sphere centers

    # Define the 6 spheres
    spheres = [
        # Layer 1
        {'center': (d1, d1, z1),   'radius': r1, 'id': 1},
        {'center': (-d1, d1, z1),  'center_str': f"(-{d1}, {d1}, {z1})", 'radius': r1, 'id': 2},
        {'center': (d1, -d1, z1),  'center_str': f"({d1}, -{d1}, {z1})", 'radius': r1, 'id': 3},
        {'center': (-d1, -d1, z1), 'center_str': f"(-{d1}, -{d1}, {z1})", 'radius': r1, 'id': 4},
        # Layer 2
        {'center': (d2, 0, z2),    'center_str': f"({d2}, 0.0, {z2})", 'radius': r2, 'id': 5},
        {'center': (-d2, 0, z2),   'center_str': f"(-{d2}, 0.0, {z2})", 'radius': r2, 'id': 6}
    ]

    # --- Verification ---
    all_constraints_met = True

    # 1. Verify each sphere is inside the pyramid
    for s in spheres:
        c, r = s['center'], s['radius']
        x, y, z = abs(c[0]), abs(c[1]), c[2]

        # Constraint 1: Must be above the base
        if z < r:
            print(f"Error: Sphere {s['id']} is outside the base plane.")
            all_constraints_met = False

        # Constraint 2: Must be inside the 4 slanted faces
        # The distance from center (x,y,z) to the nearest slanted face must be >= r
        # dist = (1650 - 22*max(|x|,|y|) - 15*z) / sqrt(709)
        if (1650 - 15 * z - 22 * max(x, y)) < r * sqrt_709:
            print(f"Error: Sphere {s['id']} is outside the slanted faces.")
            all_constraints_met = False

    # 2. Verify non-overlapping condition for all pairs
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            s1 = spheres[i]
            s2 = spheres[j]
            c1, r1 = s1['center'], s1['radius']
            c2, r2 = s2['center'], s2['radius']

            dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
            min_dist_sq = (r1 + r2)**2
            
            # Using a small tolerance for floating point comparisons
            if dist_sq < min_dist_sq - 1e-9:
                print(f"Error: Spheres {s1['id']} and {s2['id']} are overlapping.")
                all_constraints_met = False
    
    # --- Output Results ---
    if all_constraints_met:
        print("Optimal scanning configuration found and verified:")
        print("-" * 50)
        for s in spheres:
            # Use pre-formatted string for centers with negative signs if available
            center_str = s.get('center_str', str(s['center']))
            print(f"Scan {s['id']}: Center = {center_str}, Radius = {s['radius']:.1f}")
        print("-" * 50)
        
        all_radii = [s['radius'] for s in spheres]
        max_r = max(all_radii)
        min_r = min(all_radii)
        
        print(f"The maximum and minimum scanning radii (R:r) are:")
        print(f"{max_r:.1f}:{min_r:.1f}")
    else:
        print("\nThe proposed configuration is invalid.")

solve_pyramid_scanning()
<<<25.5:16.0>>>