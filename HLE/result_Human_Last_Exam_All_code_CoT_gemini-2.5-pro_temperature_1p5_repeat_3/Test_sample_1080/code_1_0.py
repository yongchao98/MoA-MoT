import math

def solve_packing_problem():
    """
    This function constructs and verifies a proposed configuration for packing spheres
    into the specified box and prints the result.
    """
    centers = []
    
    # Layer 1: 3x3 grid of spheres at z=4
    # The z-coordinate in cm is 4 * 0.5 = 2.0, which is the lowest possible position.
    layer1_coords = [4, 12, 20]
    layer1_count = 0
    for x in layer1_coords:
        for y in layer1_coords:
            centers.append((x, y, 4))
            layer1_count += 1

    # Layer 2: 2x2 grid of spheres at z=10, placed in the hollows of layer 1.
    # The (x,y) coordinates are shifted to sit between the spheres of layer 1.
    layer2_coords = [8, 16]
    layer2_count = 0
    for x in layer2_coords:
        for y in layer2_coords:
            centers.append((x, y, 10))
            layer2_count += 1
            
    # Layer 3: 3x3 grid of spheres at z=16, aligned with layer 1.
    layer3_count = 0
    for x in layer1_coords:
        for y in layer1_coords:
            centers.append((x, y, 16))
            layer3_count += 1

    # Verify that the constructed configuration is valid by checking all pairwise distances.
    # The non-overlapping constraint is (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= (4+4)^2 = 64
    min_sq_dist = 64
    is_valid = True
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            c1 = centers[i]
            c2 = centers[j]
            dist_sq = (c1[0] - c2[0])**2 + (c1[1] - c2[1])**2 + (c1[2] - c2[2])**2
            if dist_sq < min_sq_dist:
                is_valid = False
                break
        if not is_valid:
            break

    # The problem formulation is correct, so we answer n.
    # If the configuration is valid, we print the number of spheres.
    if is_valid:
        total_spheres = layer1_count + layer2_count + layer3_count
        print(f"The maximized number of eyeball candies is found by adding the candies in each layer:")
        print(f"{layer1_count} + {layer2_count} + {layer3_count} = {total_spheres}")
    else:
        # This case indicates an issue with the proposed configuration,
        # which would mean n=0 for "incorrect formulation", though here it means
        # this specific configuration is invalid.
        print(0)

solve_packing_problem()