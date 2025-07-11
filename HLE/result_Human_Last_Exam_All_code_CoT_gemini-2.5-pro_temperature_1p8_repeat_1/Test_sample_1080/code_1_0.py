def solve_packing_problem():
    """
    This function calculates the maximum number of spheres based on a dense,
    structured packing arrangement and verifies its validity.
    """
    
    # The non-overlapping constraint is that the squared distance 
    # between sphere centers must be >= 64.
    min_dist_sq = 64

    # This list will hold the (x, y, z) integer coordinates for the center of each sphere.
    centers = []

    # Layer 1: A 3x3 grid of 9 spheres at z=4.
    # The x and y coordinates are spaced by 8 units (4cm), which is the minimum distance.
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 4))

    # Layer 2: A 2x2 grid of 4 spheres at z=11, placed in the hollows of Layer 1.
    # These centers are offset in x and y to fit between the spheres below.
    for x in [8, 16]:
        for y in [8, 16]:
            centers.append((x, y, 11))

    # Layer 3: A 3x3 grid of 9 spheres at z=18, aligned with Layer 1.
    for x in [4, 12, 20]:
        for y in [4, 12, 20]:
            centers.append((x, y, 18))

    # Verify that the constructed packing is valid.
    # Check every pair of spheres for the non-overlapping constraint.
    is_valid = True
    for i in range(len(centers)):
        for j in range(i + 1, len(centers)):
            p1 = centers[i]
            p2 = centers[j]
            dist_sq = (p1[0] - p2[0])**2 + (p1[1] - p2[1])**2 + (p1[2] - p2[2])**2
            
            if dist_sq < min_dist_sq:
                # If any pair is too close, the configuration is invalid.
                is_valid = False
                break
        if not is_valid:
            break

    # If the configuration is valid, print the total number of spheres.
    # Otherwise, print 0, indicating a failure in the proposed configuration.
    if is_valid:
        # The final answer n is the total number of spheres placed.
        n = len(centers)
        print(n)
    else:
        # This case should not be reached with the correct configuration.
        print(0)

# Execute the function to find and print the maximized value n.
solve_packing_problem()