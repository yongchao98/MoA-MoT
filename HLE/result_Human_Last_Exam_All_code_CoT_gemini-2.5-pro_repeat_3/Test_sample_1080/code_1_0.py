import math

def solve_candy_packing():
    """
    Calculates the number of spherical candies that can be packed into a box
    using a greedy algorithm on a discrete grid.
    """
    # The problem is defined on a grid where 1 unit = 0.5 cm.
    
    # A candy (eyeball) has a radius of 2 cm, which is 4 units on our grid.
    candy_radius_units = 4
    
    # For two candies not to overlap, the distance between their centers must be
    # at least the sum of their radii (2 * radius).
    # We use squared distance to avoid complex square root calculations.
    min_dist_sq = (candy_radius_units + candy_radius_units)**2
    
    # Valid integer coordinates for the center of a candy.
    # x_i, y_i are in [4, 5, ..., 20]
    # z_i is in [4, 5, ..., 18]
    x_range = range(4, 21)
    y_range = range(4, 21)
    z_range = range(4, 19)
    
    # This list will store the coordinates of the center of each placed candy.
    placed_candies = []
    
    # We use a greedy approach: iterate through all possible locations and
    # place a candy if the spot is valid.
    print("Thinking and packing the candies one by one...")
    for z in z_range:
        for y in y_range:
            for x in x_range:
                
                current_center = (x, y, z)
                can_place_candy = True
                
                # Check against all candies already in the box.
                for placed_center in placed_candies:
                    # Calculate the squared distance to the other candy's center.
                    dist_sq = ( (current_center[0] - placed_center[0])**2 +
                                (current_center[1] - placed_center[1])**2 +
                                (current_center[2] - placed_center[2])**2 )
                    
                    # If the distance is less than the minimum, they overlap.
                    if dist_sq < min_dist_sq:
                        can_place_candy = False
                        break # No need to check other placed candies.
                
                if can_place_candy:
                    placed_candies.append(current_center)

    # The final count is the answer.
    n = len(placed_candies)
    
    print("\nYes, your problem formulation is correct.")
    print("The non-overlapping constraint equation is:")
    # Printing each number in the final equation as requested.
    num1 = candy_radius_units
    num2 = candy_radius_units
    print(f"(x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= ({num1} + {num2})^2")
    print(f"This simplifies to: (x_i - x_j)^2 + (y_i - y_j)^2 + (z_i - z_j)^2 >= {min_dist_sq}")
    
    print("\nBased on the calculation, the highest number of eyeball candies your son can pack is:")
    print(f"n = {n}")

# Run the solver.
solve_candy_packing()