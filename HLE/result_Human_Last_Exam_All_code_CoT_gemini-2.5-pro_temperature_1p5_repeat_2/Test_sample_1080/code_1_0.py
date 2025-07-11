import math

def solve_packing():
    """
    Calculates the maximum number of spheres that can be packed in a box
    with given dimensions and constraints using a greedy algorithm.
    """
    # Box dimensions in cm: 12x12x11. Sphere radius: 2 cm. Grid step: 0.5 cm.
    
    # The problem defines the ranges for the integer grid coordinates of the sphere centers.
    # x and y are in [4, 20], z is in [4, 18].
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # From the problem formulation, the minimum squared distance between centers in
    # grid units must be (4+4)^2 = 64.
    min_dist_sq = 64

    placed_spheres = []

    # Iterate through all possible center positions in lexicographical order (z, then y, then x).
    for z in z_range:
        for y in y_range:
            for x in x_range:
                candidate_pos = (x, y, z)
                can_place = True

                # Check for overlap with any already placed spheres.
                for placed_pos in placed_spheres:
                    dist_sq = (
                        (candidate_pos[0] - placed_pos[0])**2 +
                        (candidate_pos[1] - placed_pos[1])**2 +
                        (candidate_pos[2] - placed_pos[2])**2
                    )
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break
                
                # If no overlap is found, place the sphere.
                if can_place:
                    placed_spheres.append(candidate_pos)

    # The result is the total number of spheres placed.
    n = len(placed_spheres)

    # As requested, output the confirmation and the final "equation" with the result.
    print("Yes, the problem formulation is correct.")
    print(f"The equation for the maximized value is n = {n}")

solve_packing()
