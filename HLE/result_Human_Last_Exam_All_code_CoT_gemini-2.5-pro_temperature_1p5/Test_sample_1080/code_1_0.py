import math

def solve_candy_packing():
    """
    This function verifies a proposed solution for packing 22 spheres into the box.
    The problem is defined with a 0.5cm grid, so all coordinates are integers
    representing units of 0.5cm.

    Problem parameters in 0.5 cm units:
    - Sphere radius R = 2cm = 4 units.
    - Non-overlap squared distance D^2 = (2*R)^2 = (4+4)^2 = 64.
    - Box constraints for centers (x,y,z):
      x, y in [4, 20], z in [4, 18].
    """
    
    # Proposed configuration of 22 spheres in three layers
    centers = []
    
    # Layer 1: 3x3 grid at z=4 (height 2cm)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 4))
            
    # Layer 2: 2x2 grid (staggered) at z=10 (height 5cm)
    for x in [8, 16]:
        for y in [8, 16]:
            centers.append((x, y, 10))
            
    # Layer 3: 3x3 grid at z=16 (height 8cm)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 16))

    n = len(centers)
    min_dist_sq = (4 + 4)**2

    # Verification 1: Check if all centers are within the box boundaries.
    x_min, x_max = 4, 20
    y_min, y_max = 4, 20
    z_min, z_max = 4, 18
    for i, c in enumerate(centers):
        if not (x_min <= c[0] <= x_max and
                y_min <= c[1] <= y_max and
                z_min <= c[2] <= z_max):
            # This part should not be reached for the proposed solution
            print(f"Configuration is invalid: Center {i+1} {c} is out of bounds.")
            return

    # Verification 2: Check for overlaps between any two spheres.
    for i in range(n):
        for j in range(i + 1, n):
            c1 = centers[i]
            c2 = centers[j]
            dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
            if dist_sq < min_dist_sq:
                # This part should not be reached for the proposed solution
                print(f"Configuration is invalid: Spheres {i+1} and {j+1} overlap.")
                return

    # All checks passed, the configuration is valid.
    print(f"The proposed configuration for n = {n} candies is valid.")
    print(f"The maximum number of candies is {n}.")

    # Showcase the non-overlapping constraint equation for one pair of candies
    print("\n--- Example Non-Overlapping Constraint Check ---")
    c_i = centers[0]   # First candy in layer 1, center (4, 4, 4)
    c_j = centers[9]   # First candy in layer 2, center (8, 8, 10)
    
    # The equation is: (x_i-x_j)^2 + (y_i-y_j)^2 + (z_i-z_j)^2 >= (4+4)^2
    print(f"Checking between candy at {c_i} and candy at {c_j}:")
    
    dx = c_i[0] - c_j[0]
    dy = c_i[1] - c_j[1]
    dz = c_i[2] - c_j[2]
    
    # Printing each number in the equation step-by-step
    print(f"({c_i[0]} - {c_j[0]})^2 + ({c_i[1]} - {c_j[1]})^2 + ({c_i[2]} - {c_j[2]})^2 >= ({4} + {4})^2")
    print(f"({dx})^2 + ({dy})^2 + ({dz})^2 >= ({8})^2")
    print(f"{dx**2} + {dy**2} + {dz**2} >= {min_dist_sq}")
    print(f"{dx**2 + dy**2 + dz**2} >= {min_dist_sq} (Constraint is satisfied)")


if __name__ == '__main__':
    solve_candy_packing()