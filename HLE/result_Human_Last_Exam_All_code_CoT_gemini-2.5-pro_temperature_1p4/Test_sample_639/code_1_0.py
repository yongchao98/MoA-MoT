import math

def calculate_hard_sphere_distances():
    """
    Calculates and prints the possible normalized distances (r) between two hard
    spheres in a 2D arrangement for r <= 3.

    The distance r is normalized by the sphere diameter d.
    We assume d=1 for simplicity in calculations.
    """
    print("This script calculates the set of possible distances between two hard spheres on a 2D plane for r <= 3.")
    print("The distance 'r' is normalized by the sphere diameter 'd'.\n")

    # The normalized diameter is 1.0
    d = 1.0
    
    # Use a set to store unique distances to avoid duplicates
    distances = set()

    # --- 1. Direct Contact ---
    print("Distance 1: Two spheres in direct contact.")
    r1 = d
    distances.add(r1)
    print(f"   Equation: r = d")
    print(f"   Calculation: r = {d:.2f}")
    print(f"   Result: {r1:.2f}\n")

    # --- 2. Rhombus Configuration ---
    print("Distance 2: Long diagonal of a rhombus formed by four touching spheres.")
    val_sqrt_3 = 3.0
    r2 = math.sqrt(val_sqrt_3) * d
    distances.add(r2)
    print(f"   Equation: r = sqrt(3) * d")
    print(f"   Calculation: r = sqrt({val_sqrt_3:.2f}) * {d:.2f}")
    print(f"   Result: {r2:.2f}\n")

    # --- 3. Linear Chain of 3 ---
    print("Distance 3: Two outer spheres in a linear chain of three.")
    val_2 = 2.0
    r3 = val_2 * d
    distances.add(r3)
    print(f"   Equation: r = 2 * d")
    print(f"   Calculation: r = {val_2:.2f} * {d:.2f}")
    print(f"   Result: {r3:.2f}\n")

    # --- 4. Zig-Zag Chain of 4 ---
    print("Distance 4: First and fourth spheres in a zig-zag chain.")
    # This distance is derived from coordinate geometry: sqrt((2.5d)^2 + (d*sqrt(3)/2)^2)
    val_sqrt_7 = 7.0
    r4 = math.sqrt(val_sqrt_7) * d
    distances.add(r4)
    print(f"   Equation: r = sqrt(7) * d")
    print(f"   Calculation: r = sqrt({val_sqrt_7:.2f}) * {d:.2f}")
    print(f"   Result: {r4:.2f}\n")

    # --- 5. Linear Chain of 4 ---
    print("Distance 5: Two outer spheres in a linear chain of four.")
    val_3 = 3.0
    r5 = val_3 * d
    distances.add(r5)
    print(f"   Equation: r = 3 * d")
    print(f"   Calculation: r = {val_3:.2f} * {d:.2f}")
    print(f"   Result: {r5:.2f}\n")

    # --- Final Result ---
    sorted_distances = sorted(list(distances))
    print("-----------------------------------------------------------------")
    print("The final sorted set of unique distances r <= 3 is:")
    
    # Format the final list as requested
    final_set = [f"{dist:.2f}" for dist in sorted_distances]
    print(final_set)

if __name__ == '__main__':
    calculate_hard_sphere_distances()

<<<['1.00', '1.73', '2.00', '2.65', '3.00']>>>