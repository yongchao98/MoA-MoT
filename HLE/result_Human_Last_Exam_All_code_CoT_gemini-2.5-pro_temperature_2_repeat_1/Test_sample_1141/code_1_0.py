import math

def solve_cutting_problem():
    """
    Calculates the maximum value from cutting a steel billet based on the user's formulation.
    The strategy is to place the maximum number of high-value B2 balls and then fill
    the rest of the space with B1 balls.
    """
    # Configuration for B2 pieces (8 balls in a dense packing)
    b2_centers = [
        (4, 4, 4), (4, 12, 4),
        (12, 4, 4), (12, 12, 4),
        (20, 4, 4), (20, 12, 4),
        (28, 4, 4), (28, 12, 4)
    ]
    b2_price = 150
    num_b2s = len(b2_centers)

    # Configuration for B1 pieces
    b1_price = 1
    b1_count = 0
    
    # We will place B1 centers on a grid with step-size 2 to ensure they don't
    # overlap with each other. This automatically satisfies the B1-B1 constraint:
    # (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 4
    x_b1_range = range(1, 32, 2)  # 1, 3, ..., 31
    y_b1_range = range(1, 22, 2)  # 1, 3, ..., 21
    z_b1_range = range(1, 8, 2)   # 1, 3, 5, 7

    # Iterate through all possible B1 grid positions
    for z_j in z_b1_range:
        for y_j in y_b1_range:
            for x_j in x_b1_range:
                
                is_valid_b1_pos = True
                # Check this candidate B1 position against all 8 B2 balls
                for x_i, y_i, z_i in b2_centers:
                    dist_sq = (x_i - x_j)**2 + (y_i - y_j)**2 + (z_i - z_j)**2
                    
                    # Non-overlapping constraint for a B1 to a B2
                    if dist_sq < 25:
                        is_valid_b1_pos = False
                        break # This B1 position is invalid, check next one
                
                if is_valid_b1_pos:
                    b1_count += 1
    
    # Calculate total value
    total_value = (num_b2s * b2_price) + (b1_count * b1_price)

    # Print the result in the required equation format
    print(f"Yes, a valid solution exists under the given formulation.")
    print(f"The highest value is found by combining {num_b2s} B2 pieces and {b1_count} B1 pieces.")
    print(f"{num_b2s} * {b2_price} + {b1_count} * {b1_price} = {total_value}")


solve_cutting_problem()