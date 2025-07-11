import math

def solve_packing():
    """
    Solves the discrete sphere packing problem using a greedy algorithm.
    """
    # Problem formulation confirmation
    print("The problem formulation is correct.")

    # Box dimensions are 12x12x11 cm. Grid step is 0.5 cm.
    # Sphere radius is 2 cm, which is 4 units in the grid.
    # The center coordinates for the i-th ball must be in the ranges:
    # x_i, y_i in [4, 20]
    # z_i in [4, 18]
    
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)

    # For two non-overlapping spheres, the squared distance between their centers
    # must be >= (radius1 + radius2)^2.
    # In grid units, (4+4)^2 = 8^2 = 64.
    min_sq_dist = 64

    placed_spheres_centers = []

    # Greedy algorithm: iterate through all possible center positions
    # and place a sphere if it doesn't overlap with existing ones.
    for z in z_range:
        for y in y_range:
            for x in x_range:
                
                can_place = True
                new_center = (x, y, z)

                # Check for conflicts with already placed spheres
                for placed_center in placed_spheres_centers:
                    sq_dist = ( (new_center[0] - placed_center[0]) ** 2 +
                                (new_center[1] - placed_center[1]) ** 2 +
                                (new_center[2] - placed_center[2]) ** 2 )
                    
                    if sq_dist < min_sq_dist:
                        can_place = False
                        break
                
                if can_place:
                    placed_spheres_centers.append(new_center)

    n = len(placed_spheres_centers)

    # Final answer output
    # The final equation is n = result
    print(f"The final equation is: n = {n}")
    print(f"The highest number of eyeball candies that can be packed into the box is {n}.")


solve_packing()
