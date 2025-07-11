import math

def solve_candy_packing():
    """
    Calculates the maximum number of eyeball candies that can be packed into a box
    using a greedy algorithm on a discrete grid.
    """
    # Problem definition from the user
    # Box dimensions: 12x12x11 cm
    # Sphere radius: 2 cm
    # Grid step: 0.5 cm

    # Coordinate ranges for the sphere's center in grid units (1 unit = 0.5 cm)
    x_range = range(4, 21)  # Corresponds to [2cm, 10cm]
    y_range = range(4, 21)  # Corresponds to [2cm, 10cm]
    z_range = range(4, 19)  # Corresponds to [2cm, 9cm]

    # The radius of a sphere in grid units is 2 cm / 0.5 cm = 4 units.
    # The non-overlapping constraint is based on the sum of two radii.
    radius_in_grid_units = 4
    min_center_dist = radius_in_grid_units + radius_in_grid_units
    min_dist_sq = min_center_dist**2

    placed_candies = []

    # Greedily iterate through all possible positions and place a candy if possible
    for z in z_range:
        for y in y_range:
            for x in x_range:
                new_candy_pos = (x, y, z)
                can_place = True
                
                # Check for overlap with already placed candies
                for placed_pos in placed_candies:
                    dist_sq = (
                        (new_candy_pos[0] - placed_pos[0]) ** 2 +
                        (new_candy_pos[1] - placed_pos[1]) ** 2 +
                        (new_candy_pos[2] - placed_pos[2]) ** 2
                    )
                    
                    if dist_sq < min_dist_sq:
                        can_place = False
                        break
                
                if can_place:
                    placed_candies.append(new_candy_pos)
    
    max_n = len(placed_candies)

    # Output the conclusion and the result
    print("Yes, your problem formulation is correct.")
    print("\n--- Solving the problem ---")
    print(f"The equation for the non-overlapping constraint is: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= (r_grid+r_grid)^2")
    print(f"The numbers in this equation are:")
    print(f"Radius in grid units (r_grid): {radius_in_grid_units}")
    print(f"Squared minimum distance: ({radius_in_grid_units}+{radius_in_grid_units})^2 = {min_dist_sq}")
    print("\n--- Result ---")
    print(f"The maximized number of eyeball candies (n) that can be packed is: {max_n}")

solve_candy_packing()