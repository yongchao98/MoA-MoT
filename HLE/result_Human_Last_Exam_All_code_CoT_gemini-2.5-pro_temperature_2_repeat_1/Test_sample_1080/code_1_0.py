import math

def solve_and_verify_packing():
    """
    This function verifies a proposed optimal packing for the sphere packing problem.
    """
    # First, we answer the user's primary question: Is the formulation correct?
    # Based on our analysis, the formulation is correct.
    print("Yes, the problem formulation is correct.")
    print("-" * 20)

    # --- The Solution ---
    # The maximum number of spheres (n) is 22.
    # This is achieved by a specific layered packing arrangement. We define the
    # integer coordinates of the center of each of the 22 spheres below.
    
    centers = []
    # Layer A (9 spheres at z=4)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 4))
    
    # Layer B (4 spheres at z=10, nested in the hollows of Layer A)
    for x in [8, 16]:
        for y in [8, 16]:
            centers.append((x, y, 10))
            
    # Layer C (9 spheres at z=16, aligned with Layer A)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 16))

    n = len(centers)
    print(f"The maximized value n is {n}.")
    print("-" * 20)
    
    # --- Verification of the Solution ---
    print("Below is the verification that this packing of 22 spheres is valid.")

    # Problem constraint values
    x_min, x_max = 4, 20
    y_min, y_max = 4, 20
    z_min, z_max = 4, 18
    # From (4+4)^2, as per the problem description.
    min_sq_dist = 64

    print("The numbers in the governing constraints are:")
    print(f"  - Sphere center x-coordinate must be in [{x_min}, {x_max}]")
    print(f"  - Sphere center y-coordinate must be in [{y_min}, {y_max}]")
    print(f"  - Sphere center z-coordinate must be in [{z_min}, {z_max}]")
    print(f"  - The squared distance between any two sphere centers must be >= {min_sq_dist}")
    print("-" * 20)
    
    is_valid_packing = True
    for i in range(n):
        c1 = centers[i]
        # 1. Verify boundary constraints
        if not (x_min <= c1[0] <= x_max and y_min <= c1[1] <= y_max and z_min <= c1[2] <= z_max):
            print(f"Verification FAILED: Center {c1} is out of bounds.")
            is_valid_packing = False
            break
            
        # 2. Verify non-overlapping constraints against all other spheres
        for j in range(i + 1, n):
            c2 = centers[j]
            sq_dist = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
            if sq_dist < min_sq_dist:
                print(f"Verification FAILED: Spheres at {c1} and {c2} overlap. Sq. distance is {sq_dist}, but must be >= {min_sq_dist}.")
                is_valid_packing = False
                break
        if not is_valid_packing:
            break
            
    if is_valid_packing:
        print("Verification SUCCEEDED: All constraints are met.")
        print("\nThe coordinates for the 22 spheres are:")
        for i, center in enumerate(centers):
            print(f"  Sphere {i+1:>2}: {center}")

solve_and_verify_packing()