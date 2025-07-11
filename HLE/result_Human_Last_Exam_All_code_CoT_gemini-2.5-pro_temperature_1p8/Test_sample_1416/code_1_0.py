import math

def solve_pyramid_scanning():
    """
    Calculates the optimal radii and locations for 6 scans inside a pyramid
    based on a symmetric '4+2' configuration.
    """

    # 1. Define pyramid and scanner constants
    H = 110.0  # Pyramid height in meters
    A = 150.0  # Pyramid base side length in meters
    MIN_R = 10.0
    # Pre-calculated constant for distance formula to side faces
    K = math.sqrt(4 * H**2 + A**2)

    print("This script calculates an optimal solution for placing 6 non-overlapping spheres inside a square pyramid.")
    print("The chosen configuration is a symmetric '4+2' arrangement.\n")

    # 2. Calculate parameters for the 4 outer spheres
    print("--- Step 1: Calculating the 4 Outer Spheres ---")
    print("To maximize their size, the spheres are placed in a square pattern, as low as possible.")
    print("Their centers are assumed to be at coordinates (+/-d, +/-d, z), with z=d for compactness.")
    # For a sphere to be tangent to the pyramid base, a side wall, and its neighbors simultaneously,
    # we can set its radius 'r' to be equal to its coordinate distances 'd' and 'z'.
    # So, center is at (+/-r, +/-r, r).
    # The sphere must be inside the pyramid, so its distance to the side faces must be >= r.
    # The governing constraint is: (A*H - 2*H*r - A*r) / K >= r
    # Solving for r: r <= (A * H) / (K + 2*H + A)
    r1_max_continuous = (A * H) / (K + 2*H + A)
    # Radii must be a multiple of 0.5, so we round down.
    r1 = math.floor(r1_max_continuous / 0.5) * 0.5
    d1 = r1
    z1 = r1

    print(f"Maximum theoretical radius for these spheres: {r1_max_continuous:.2f} m.")
    print(f"Using a valid radius step, the optimal radius r1 = {r1:.1f} m.")

    outer_spheres_radii = [r1] * 4
    outer_spheres_centers = [
        (d1, d1, z1),
        (-d1, d1, z1),
        (d1, -d1, z1),
        (-d1, -d1, z1),
    ]
    print(f"Result: 4 spheres with radius {r1:.1f} m.\n")


    # 3. Calculate parameters for the 2 inner spheres
    print("--- Step 2: Calculating the 2 Inner Spheres ---")
    # These 4 spheres form a central vertical tunnel. The top of these spheres is at z = z1 + r1
    tunnel_height = z1 + r1
    print(f"The outer spheres create a central tunnel up to height z = {tunnel_height:.1f} m.")
    
    # We can fit a stack of two spheres inside this tunnel.
    # The sum of their diameters must not exceed the tunnel height: 2*r2 + 2*r3 <= tunnel_height
    # Therefore, r2 + r3 <= tunnel_height / 2
    sum_radii_inner = tunnel_height / 2
    print(f"The sum of the radii for the two inner spheres (r2 + r3) must be at most {sum_radii_inner:.1f} m.")
    print(f"To maximize their volume (proportional to r^3), we should make their radii as different as possible.")

    # We set one radius to the minimum allowed (10.0m) to maximize the other.
    r2 = MIN_R
    r3 = sum_radii_inner - r2
    
    inner_spheres_radii = [r2, r3]

    # Calculate their center positions for a vertical stack starting from z=0
    z_center_2 = r2
    z_center_3 = z_center_2 + r2 + r3

    inner_spheres_centers = [
        (0.0, 0.0, z_center_2),
        (0.0, 0.0, z_center_3),
    ]
    print(f"Result: 2 spheres with radii {r2:.1f} m and {r3:.1f} m.\n")

    # 4. Consolidate and present the final results
    print("--- Step 3: Final Configuration ---")
    all_radii = outer_spheres_radii + inner_spheres_radii
    all_centers = outer_spheres_centers + inner_spheres_centers

    print("Optimal scanning locations (x,y,z) and radii (r):")
    for i in range(len(all_radii)):
        center = all_centers[i]
        radius = all_radii[i]
        print(f"  Scan {i+1}: Center=({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}), Radius={radius:.1f}")

    max_radius = max(all_radii)
    min_radius = min(all_radii)
    
    print("\nThe question asks for the answer in the form R:r")
    print(f"Maximum radius (R) = {max_radius:.1f} m")
    print(f"Minimum radius (r) = {min_radius:.1f} m")

solve_pyramid_scanning()
<<<25.5:10.0>>>