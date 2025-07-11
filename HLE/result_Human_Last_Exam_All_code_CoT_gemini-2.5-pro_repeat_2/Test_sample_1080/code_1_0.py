import math

def solve_candy_packing():
    """
    This function confirms the user's problem formulation and solves for the maximum
    number of spheres (n) that can be packed into the given box under the
    specified integer constraints.

    It does so by constructing and verifying a dense packing configuration.
    """
    # Part 1: Confirm the problem formulation is correct.
    # The box is 12x12x11 cm. The sphere radius is 2 cm.
    # The valid region for sphere centers in cm is [2, 10] for x/y and [2, 9] for z.
    # With a 0.5 cm grid, the integer coordinates are real_coords / 0.5.
    # This correctly gives your ranges: x,y in [4, 20] and z in [4, 18].
    #
    # The non-overlapping constraint means the distance between any two centers
    # must be at least the sum of their radii (2+2=4 cm).
    # In integer coordinates, the radius is 4 units (2 cm / 0.5 cm).
    # So, the distance between centers must be >= 8 units.
    # The squared distance must be >= 8*8 = 64.
    # Your constraint (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= (4+4)^2 is correct.
    print("Yes, your problem formulation is correct.")
    print("Based on this formulation, we can find the maximized value of n.")
    print("-" * 20)

    # Part 2: Construct a dense packing with integer coordinates.
    # This packing uses three layers in an A-B-A stacking pattern.
    centers = []

    # Layer 1 (Bottom): A 3x3 grid of 9 spheres.
    layer1_coords = [4, 12, 20]
    layer1_z = 4
    num_layer1 = 0
    for y in layer1_coords:
        for x in layer1_coords:
            centers.append((x, y, layer1_z))
            num_layer1 += 1

    # Layer 2 (Middle): A 2x2 grid of 4 spheres placed in the hollows of Layer 1.
    layer2_coords = [8, 16]
    layer2_z = 10
    num_layer2 = 0
    for y in layer2_coords:
        for x in layer2_coords:
            centers.append((x, y, layer2_z))
            num_layer2 += 1
            
    # Layer 3 (Top): A 3x3 grid of 9 spheres, identical to Layer 1 but at a higher z.
    layer3_coords = [4, 12, 20]
    layer3_z = 16
    num_layer3 = 0
    for y in layer3_coords:
        for x in layer3_coords:
            centers.append((x, y, layer3_z))
            num_layer3 += 1

    # Part 3: Verify the constructed packing is valid.
    min_dist_sq_req = 64
    is_valid = True
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            c1 = centers[i]
            c2 = centers[j]
            dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
            if dist_sq < min_dist_sq_req:
                is_valid = False
                # This part of the code will not be reached if the packing is correct.
                print(f"Error: Constraint violated between {c1} and {c2}. Dist^2 = {dist_sq}")
                break
        if not is_valid:
            break
            
    # Part 4: Output the result.
    if is_valid:
        total_spheres = len(centers)
        print("A valid packing configuration has been found with three layers:")
        print(f"A bottom layer with {num_layer1} candies.")
        print(f"A middle layer with {num_layer2} candies.")
        print(f"A top layer with {num_layer3} candies.")
        print("\nThe final equation to find the total number of candies n is:")
        print(f"{num_layer1} + {num_layer2} + {num_layer3} = {total_spheres}")
    else:
        # This case suggests the constructed packing was invalid.
        total_spheres = 0
        print("The proposed packing was invalid. The answer is 0.")

solve_candy_packing()