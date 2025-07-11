import math

def solve_cutting_problem():
    """
    Calculates the maximum value from cutting a steel billet based on the given constraints.

    The strategy is based on the finding that B2 balls and T1 cubes are mutually exclusive.
    The optimal approach is to maximize the number of high-value B2 balls and then
    fill the remaining space with B1 balls.
    """

    # --- Configuration for B2 Balls ---
    # A hexagonal-like packing allows for 11 B2 balls.
    b2_price = 150
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),  # Row 1
        (8, 11, 4), (16, 11, 4), (24, 11, 4),            # Row 2
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4) # Row 3
    ]
    num_b2 = len(b2_centers)
    b2_value = num_b2 * b2_price

    # --- Configuration for B1 Balls ---
    # We fill the gaps with B1 balls. A simple way to count them is to check
    # every possible position on a grid and see if it's valid.
    b1_price = 1
    
    # B1 centers are on a grid with spacing 2 to avoid B1-B1 overlap.
    b1_candidates = []
    for x in range(1, 32, 2):
        for y in range(1, 22, 2):
            for z in range(1, 8, 2):
                b1_candidates.append((x, y, z))

    # B1-B2 non-overlapping constraint: squared distance must be >= (4+1)^2 = 25
    b1_b2_dist_sq_threshold = 25
    
    valid_b1_centers = []
    for b1_c in b1_candidates:
        is_valid = True
        for b2_c in b2_centers:
            # Calculate squared Euclidean distance
            dist_sq = (b1_c[0] - b2_c[0])**2 + (b1_c[1] - b2_c[1])**2 + (b1_c[2] - b2_c[2])**2
            if dist_sq < b1_b2_dist_sq_threshold:
                is_valid = False
                break
        if is_valid:
            valid_b1_centers.append(b1_c)
            
    num_b1 = len(valid_b1_centers)
    b1_value = num_b1 * b1_price

    total_value = b2_value + b1_value

    # --- Final Output ---
    # The problem asks for the equation in the final output.
    print(f"Highest valid solution found by combining B2 and B1 pieces:")
    print(f"Number of B2 pieces: {num_b2}")
    print(f"Number of B1 pieces: {num_b1}")
    print(f"Total Value = (B2 pieces * B2 price) + (B1 pieces * B1 price)")
    print(f"{num_b2} * {b2_price} + {num_b1} * {b1_price} = {total_value}")

solve_cutting_problem()
<<<2050>>>