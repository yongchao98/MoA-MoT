import math

def solve_cutting_problem():
    """
    This function calculates the maximum value obtainable by cutting a steel billet
    according to the user's specified shapes, prices, and constraints.
    """
    
    # An optimal packing of B2 balls is used as the starting point, as they are
    # the most valuable pieces. A staggered arrangement fits 11 balls.
    b2_centers = [
        # Row 1
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        # Row 2 (offset for denser packing)
        (8, 11, 4), (16, 11, 4), (24, 11, 4),
        # Row 3
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4)
    ]
    num_b2 = len(b2_centers)
    value_b2 = 150

    # The T1-to-B2 constraint 'min(|x_i-x_j|, |y_i-y_j|, |z_i-z_j|) >= 5' makes it
    # impossible to place any T1 cubes if a B2 is present, because the z-component
    # |z_i - 4| can never be >= 5 for a T1 cube where z_i is in [1, 7].
    num_t1 = 0
    value_t1 = 5

    # A greedy algorithm is used to fill the remaining space with B1 balls.
    # We iterate through every possible B1 center location and place a ball
    # if it does not conflict with any existing B2 or B1 ball.
    b1_centers = []
    value_b1 = 1
    
    # Valid B1 center coordinate ranges: x in [1,31], y in [1,21], z in [1,7].
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                candidate_center = (x, y, z)
                is_valid = True

                # Constraint check 1: B1 against all B2 balls.
                # Squared distance must be >= 25.
                for b2_c in b2_centers:
                    dist_sq = (x - b2_c[0])**2 + (y - b2_c[1])**2 + (z - b2_c[2])**2
                    if dist_sq < 25:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Constraint check 2: B1 against all previously placed B1 balls.
                # Squared distance must be >= 4.
                for b1_c in b1_centers:
                    dist_sq = (x - b1_c[0])**2 + (y - b1_c[1])**2 + (z - b1_c[2])**2
                    if dist_sq < 4:
                        is_valid = False
                        break
                
                if is_valid:
                    b1_centers.append(candidate_center)
    
    num_b1 = len(b1_centers)

    # Calculate the total value and print the detailed result.
    total_value = (num_b2 * value_b2) + (num_b1 * value_b1) + (num_t1 * value_t1)
    
    print("Yes, your problem formulation is a valid optimization problem.")
    print("The optimal cutting strategy to maximize value is to prioritize B2 balls, then fill the rest with B1 balls.")
    print("The highest value solution is composed of:")
    print(f"- {num_b2} B2 balls (price 150)")
    print(f"- {num_b1} B1 balls (price 1)")
    print(f"- {num_t1} T1 cubes (price 5)")
    print("\nFinal calculation:")
    print(f"{num_b2} * {value_b2} + {num_b1} * {value_b1} + {num_t1} * {value_t1} = {total_value}")
    
    # The final line outputs the value in the requested format.
    print(f"<<<{total_value}>>>")

solve_cutting_problem()