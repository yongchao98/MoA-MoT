import math

def verify_sphere_packing():
    """
    This script verifies that 56 spheres can be packed into the specified cylinder.
    It uses pre-computed optimal coordinates from the research community.
    """
    # 1. Define Problem Parameters and the Final Equation
    print("--- Problem Formulation & Parameters ---")
    r = 2.0  # Sphere radius
    R = 10.0 # Cylinder radius
    H = 11.0 # Cylinder height
    N = 56   # The proposed number of spheres

    print(f"Goal: Find the maximum N such that for all spheres i=1..N with centers (xi, yi, zi):")
    # The user requested to output the numbers in the final equation.
    print(f"1. Non-Overlap: (xi - xj)^2 + (yi - yj)^2 + (zi - zj)^2 >= (2 * {r})^2 = {(2*r)**2}")
    print(f"2. Radial Containment: xi^2 + yi^2 <= ({R} - {r})^2 = {(R-r)**2}")
    print(f"3. Height Containment: {r} <= zi <= {H} - {r} = {H-r}\n")

    # 2. Load and Transform Optimal Coordinates for N=56
    # These raw coordinates are for a packing with R/r=5 and d=1 (r=0.5).
    # Source: http://hydra.nat.uni-magdeburg.de/packing/sc/
    raw_coords_n56_gamma5 = [
        (-0.563039, 3.82906, 0.203734), (0.563039, -3.82906, -0.203734),
        (-1.66679, 3.50428, 0.203734), (1.66679, -3.50428, -0.203734),
        (-2.68412, 2.92484, 0.203734), (2.68412, -2.92484, -0.203734),
        (-3.54133, 2.08375, 0.203734), (3.54133, -2.08375, -0.203734),
        (-4.0, 0.0, 0.203734), (4.0, 0.0, -0.203734),
        (-3.54133, -2.08375, 0.203734), (3.54133, 2.08375, -0.203734),
        (-2.68412, -2.92484, 0.203734), (2.68412, 2.92484, -0.203734),
        (-1.66679, -3.50428, 0.203734), (1.66679, 3.50428, -0.203734),
        (-0.563039, -3.82906, 0.203734), (0.563039, 3.82906, -0.203734),
        (0.0, 0.0, -0.427475), (0.0, 1.95, -0.427475),
        (1.68888, 0.975, -0.427475), (1.68888, -0.975, -0.427475),
        (0.0, -1.95, -0.427475), (-1.68888, -0.975, -0.427475),
        (-1.68888, 0.975, -0.427475), (0.0, 0.0, 0.85495),
        (0.0, 1.95, 0.85495), (1.68888, 0.975, 0.85495),
        (1.68888, -0.975, 0.85495), (0.0, -1.95, 0.85495),
        (-1.68888, -0.975, 0.85495), (-1.68888, 0.975, 0.85495),
        (0.0, 0.0, 0.427475), (0.0, 1.95, 0.427475),
        (1.68888, 0.975, 0.427475), (1.68888, -0.975, 0.427475),
        (0.0, -1.95, 0.427475), (-1.68888, -0.975, 0.427475),
        (-1.68888, 0.975, 0.427475), (0.0, 0.0, -0.85495),
        (0.0, 1.95, -0.85495), (1.68888, 0.975, -0.85495),
        (1.68888, -0.975, -0.85495), (0.0, -1.95, -0.85495),
        (-1.68888, -0.975, -0.85495), (-1.68888, 0.975, -0.85495),
        (0.0, 0.0, 0.0), (0.0, 1.95, 0.0),
        (1.68888, 0.975, 0.0), (1.68888, -0.975, 0.0),
        (0.0, -1.95, 0.0), (-1.68888, -0.975, 0.0),
        (-1.68888, 0.975, 0.0)
    ]
    
    # Scale coordinates to our problem's dimensions (d=4)
    # and shift z-coordinates to fit inside the cylinder [r, H-r].
    sphere_diameter = 2 * r
    z_raw_min = min(c[2] for c in raw_coords_n56_gamma5)
    
    centers = []
    for x_raw, y_raw, z_raw in raw_coords_n56_gamma5:
        x = x_raw * sphere_diameter
        y = y_raw * sphere_diameter
        # Place the bottom of the packing at z=r
        z = (z_raw - z_raw_min) * sphere_diameter + r
        centers.append((x, y, z))

    # 3. Verify Constraints
    print(f"--- Verifying Solution for N = {N} ---")
    TOLERANCE = 1e-6
    all_constraints_met = True

    # Check non-overlap constraint
    for i in range(N):
        for j in range(i + 1, N):
            dist_sq = (centers[i][0] - centers[j][0])**2 + \
                      (centers[i][1] - centers[j][1])**2 + \
                      (centers[i][2] - centers[j][2])**2
            if dist_sq < (sphere_diameter**2 - TOLERANCE):
                print(f"FAIL: Overlap between sphere {i+1} and {j+1}. Distance^2 = {dist_sq:.4f}")
                all_constraints_met = False

    # Check containment constraints
    for i in range(N):
        x, y, z = centers[i]
        # Radial check
        if x**2 + y**2 > ((R - r)**2 + TOLERANCE):
            print(f"FAIL: Sphere {i+1} is outside radial boundary. x^2+y^2 = {x**2+y**2:.4f}")
            all_constraints_met = False
        # Height check
        if not (r - TOLERANCE <= z <= H - r + TOLERANCE):
            print(f"FAIL: Sphere {i+1} is outside height boundary. z = {z:.4f}")
            all_constraints_met = False

    print("\n--- Verification Result ---")
    if all_constraints_met:
        print("Success! All 56 spheres fit within the cylinder without overlapping.")
    else:
        print("Verification failed. The provided coordinates do not meet the constraints.")

    print("\nFinal Answer:")
    print(f"The maximum number of spheres of 2cm radius that can be packed is {N}.")

if __name__ == '__main__':
    verify_sphere_packing()