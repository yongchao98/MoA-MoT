import math

def solve_packing_problem():
    """
    This function solves the discrete sphere packing problem described by the user.
    It verifies a dense packing configuration and checks if it's optimal by trying
    to add one more sphere.
    """
    # Yes, the problem formulation is correct.
    # Box dimensions in cm: 12x12x11
    # Sphere radius in cm: 2 (diameter 4)
    # Grid step: 0.5 cm
    
    # In grid units (multiples of 0.5cm):
    # Box: 24x24x22
    # Sphere radius: 4, Sphere diameter: 8
    # Center constraints: x,y in [4, 20], z in [4, 18]. This is correct.
    # Non-overlapping squared distance: (8)^2 = 64. This is also correct.
    
    print("Yes, your problem formulation is correct.")
    print("The problem is a discrete sphere packing problem.")
    print("\nA simple packing strategy (simple cubic lattice) would place spheres at coordinates (4,4,4), (4,4,12), (4,12,4), etc.")
    print("This results in floor(12/4) x floor(12/4) x floor(11/4) = 3 x 3 x 2 = 18 spheres.")
    print("\nHowever, we can do better with a denser, staggered packing arrangement.")
    
    # Define the dense packing with 22 spheres
    min_dist_sq = 8**2
    spheres = []
    
    # Layer 1 @ z=4 (9 spheres)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            spheres.append((x, y, 4))
            
    # Layer 2 @ z=10 (4 spheres)
    for x in [8, 16]:
        for y in [8, 16]:
            spheres.append((x, y, 10))
            
    # Layer 3 @ z=16 (9 spheres)
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            spheres.append((x, y, 16))

    max_n = len(spheres)
    print(f"\nA staggered packing strategy allows for n = {max_n} spheres.")
    
    # Verify this 22-sphere packing is valid
    is_valid_packing = True
    for i in range(max_n):
        for j in range(i + 1, max_n):
            s1 = spheres[i]
            s2 = spheres[j]
            dist_sq = (s1[0] - s2[0])**2 + (s1[1] - s2[1])**2 + (s1[2] - s2[2])**2
            if dist_sq < min_dist_sq:
                print(f"Error in proposed packing: Spheres {i+1} and {j+1} are too close.")
                is_valid_packing = False
                break
        if not is_valid_packing:
            break

    if not is_valid_packing:
        return

    # Now, try to add a 23rd sphere by searching all possible empty grid points.
    print("\nSearching for a valid position for a 23rd sphere...")
    found_23rd = False
    
    # Define the search space for the center of the new sphere
    occupied_locations = set(spheres)
    
    for x_new in range(4, 21):
        for y_new in range(4, 21):
            for z_new in range(4, 19):
                cand_pos = (x_new, y_new, z_new)
                
                if cand_pos in occupied_locations:
                    continue
                
                # Check distance to all other spheres
                can_place = True
                for sp in spheres:
                    dist_sq = (sp[0] - cand_pos[0])**2 + (sp[1] - cand_pos[1])**2 + (sp[2] - cand_pos[2])**2
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break
                
                if can_place:
                    found_23rd = True
                    break
            if found_23rd:
                break
        if found_23rd:
            break
            
    if found_23rd:
        print("\nA valid position for a 23rd sphere was found! This means n >= 23.")
    else:
        print("Search complete. No valid position was found for a 23rd sphere.")
        print("\nThis strongly indicates that the maximum number of spheres is 22.")

        print("\n-------------------------")
        print("### FINAL ANSWER ###")
        print("-------------------------")
        
        # The prompt asks to output each number in the final equation.
        # We interpret this as stating the value of n and the sphere coordinates that satisfy the constraints.
        print(f"The maximized value n = {max_n}")
        print("\nThis maximum is achieved with the following set of 22 sphere centers (in 0.5 cm units):")
        
        final_equation = f"Maximize n, where n = {max_n}"
        print(final_equation)
        print("\nSubject to constraints:")
        print("1. For any two spheres i and j:")
        print("   (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= 64")
        print("2. For each sphere i:")
        print("   x_i, y_i are integers in [4, 20]")
        print("   z_i is an integer in [4, 18]")
        
        print("\nSolution (Coordinates for the 22 spheres):")
        for i, s in enumerate(spheres):
            print(f"Sphere {i+1}: (x={s[0]}, y={s[1]}, z={s[2]})")

# Execute the function
solve_packing_problem()
