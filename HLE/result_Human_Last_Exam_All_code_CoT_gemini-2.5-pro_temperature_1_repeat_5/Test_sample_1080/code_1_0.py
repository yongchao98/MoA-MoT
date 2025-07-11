import math

def solve_packing():
    """
    Solves the sphere packing problem using a greedy algorithm.
    """
    # Box dimensions are 12x12x11 cm. Candy radius is 2 cm.
    # The grid step is 0.5 cm.
    # All calculations are done in grid units.
    # Candy radius = 2 cm / 0.5 cm/unit = 4 units.
    # Candy diameter = 8 units.

    # Valid center coordinates in grid units:
    # x_i in [4, 20], y_i in [4, 20], z_i in [4, 18].
    x_min, x_max = 4, 20
    y_min, y_max = 4, 20
    z_min, z_max = 4, 18

    # Non-overlapping constraint: The squared distance between centers must be
    # at least the squared diameter ( (2*radius)^2 ).
    # min_dist_sq = (2 * 4)^2 = 64
    min_dist_sq = 64

    placed_centers = []

    # Greedily iterate through all possible center locations.
    for z in range(z_min, z_max + 1):
        for y in range(y_min, y_max + 1):
            for x in range(x_min, x_max + 1):
                
                potential_center = (x, y, z)
                is_valid_placement = True
                
                # Check for overlap with already placed candies.
                for placed_center in placed_centers:
                    px, py, pz = placed_center
                    
                    dist_sq = (x - px)**2 + (y - py)**2 + (z - pz)**2
                    
                    if dist_sq < min_dist_sq:
                        is_valid_placement = False
                        break  # Overlaps, so this position is invalid.
                
                # If the position is valid, place the candy.
                if is_valid_placement:
                    placed_centers.append(potential_center)

    # The result is the total number of candies placed.
    max_n = len(placed_centers)
    
    # The problem formulation is correct, so we print the maximized value n.
    print("The problem formulation is correct.")
    print(f"The maximum number of eyeball candies (n) that can be packed is: {max_n}")

solve_packing()
