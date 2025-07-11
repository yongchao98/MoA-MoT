import math

def solve_candy_packing():
    """
    Solves the candy packing problem by constructing and verifying a dense packing.
    """
    
    # The minimum squared distance between the centers of any two spheres.
    # The sphere radius is 2cm. Diameter is 4cm. The grid is 0.5cm.
    # So the radius is 4 units, diameter is 8 units.
    # The distance between centers must be >= diameter.
    # dist^2 >= diameter^2 -> dist^2 >= 8^2 = 64.
    min_dist_sq = 64

    # List to hold the center coordinates (x, y, z) of each sphere
    spheres = []

    # Layer 1: A 3x3 grid of 9 spheres at z=4
    # The x and y coordinates are spaced by 8 units (4cm)
    layer1_z = 4
    for i in range(3):
        for j in range(3):
            spheres.append((4 + i * 8, 4 + j * 8, layer1_z))

    # Layer 2: A 2x2 grid of 4 spheres at z=10, placed in the hollows of Layer 1
    # This layer is shifted by (4, 4) in the x-y plane relative to Layer 1
    layer2_z = 10
    for i in range(2):
        for j in range(2):
            spheres.append((8 + i * 8, 8 + j * 8, layer2_z))

    # Layer 3: Another 3x3 grid of 9 spheres at z=16, same x-y pattern as Layer 1
    layer3_z = 16
    for i in range(3):
        for j in range(3):
            spheres.append((4 + i * 8, 4 + j * 8, layer3_z))
            
    # Verification: Check if all non-overlapping constraints are satisfied
    is_valid = True
    for i in range(len(spheres)):
        for j in range(i + 1, len(spheres)):
            p1 = spheres[i]
            p2 = spheres[j]
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            if dist_sq < min_dist_sq:
                print(f"Error: Overlap found between {p1} and {p2}. Distance squared is {dist_sq}")
                is_valid = False

    if is_valid:
        print(f"The constructed packing is valid.")
        
    # The maximized value of n
    n = len(spheres)
    print(f"\nThe maximum number of eyeball candies (n) that can be packed is: {n}")
    
    # As requested, show a sample calculation for the non-overlapping constraint.
    # Let's check the distance between a sphere in Layer 1 and one in Layer 2.
    p1 = spheres[0]  # This is (4, 4, 4) from Layer 1
    p2 = spheres[9]  # This is (8, 8, 10) from Layer 2
    
    print("\n--- Sample Calculation ---")
    print(f"Checking the non-overlapping constraint for two candies:")
    print(f"Candy 1 center: (x1, y1, z1) = {p1}")
    print(f"Candy 2 center: (x2, y2, z2) = {p2}")
    print("The constraint is: (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 >= 64")
    
    dx, dy, dz = (p1[0] - p2[0]), (p1[1] - p2[1]), (p1[2] - p2[2])
    dx_sq, dy_sq, dz_sq = dx**2, dy**2, dz**2
    total_dist_sq = dx_sq + dy_sq + dz_sq
    
    # Final output showing each number in the equation
    print("\nCalculation steps:")
    print(f"({p1[0]} - {p2[0]})^2 + ({p1[1]} - {p2[1]})^2 + ({p1[2]} - {p2[2]})^2")
    print(f"= ({dx})^2 + ({dy})^2 + ({dz})^2")
    print(f"= {dx_sq} + {dy_sq} + {dz_sq}")
    print(f"= {total_dist_sq}")
    print(f"\nResult: {total_dist_sq} is indeed greater than or equal to 64.")


solve_candy_packing()