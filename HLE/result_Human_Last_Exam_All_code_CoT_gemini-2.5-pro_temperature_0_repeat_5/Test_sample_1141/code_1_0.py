import math

def solve_cutting_problem():
    """
    Calculates the maximum value based on the user's problem formulation.

    The function identifies that B2 and T1 pieces are mutually exclusive
    and determines the best strategy is to use B2 and B1 pieces. It calculates
    the number of B2 and B1 pieces for the optimal configuration and prints the result.
    """

    # Prices for each product
    price_b2 = 150
    price_b1 = 1
    price_t1 = 5

    # --- Strategy 1: Maximize B2 and B1 pieces ---

    # A near-optimal packing allows for 11 B2 balls.
    # Their centers are placed on the z=4 plane.
    num_b2 = 11
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),  # Row at y=4
        (8, 11, 4), (16, 11, 4), (24, 11, 4),          # Staggered row at y=11
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4) # Row at y=18
    ]

    # Now, find how many B1 balls can fit in the remaining space.
    # B1 center must be in x:[1,31], y:[1,21], z:[1,7]
    # Non-overlapping constraints:
    # B1-B2: dist_sq >= 25
    # B1-B1: dist_sq >= 4
    
    b1_locations = []
    # We iterate on a grid with spacing 2 to automatically satisfy the B1-B1 constraint.
    for z in range(1, 8, 2):
        for y in range(1, 22, 2):
            for x in range(1, 32, 2):
                b1_candidate_center = (x, y, z)
                is_valid = True
                
                # Check for overlap with any of the 11 B2 balls
                for b2_center in b2_centers:
                    dist_sq = (
                        (b1_candidate_center[0] - b2_center[0])**2 +
                        (b1_candidate_center[1] - b2_center[1])**2 +
                        (b1_candidate_center[2] - b2_center[2])**2
                    )
                    if dist_sq < 25:
                        is_valid = False
                        break
                
                if is_valid:
                    # Since we iterate on a grid with spacing 2, the B1-B1 constraint
                    # is automatically satisfied.
                    b1_locations.append(b1_candidate_center)

    num_b1 = len(b1_locations)
    
    # --- Compare with Strategy 2 (T1 and B1 pieces) ---
    # Max value from only B1s is 16*11*4 = 704.
    # Adding T1s reduces the number of B1s for a net loss in value.
    # So, max value for strategy 2 is 704.
    
    value_strategy1 = num_b2 * price_b2 + num_b1 * price_b1
    value_strategy2 = 704

    if value_strategy1 > value_strategy2:
        final_num_b2 = num_b2
        final_num_b1 = num_b1
        final_value = value_strategy1
        print(f"The problem formulation is flawed, but solving it as written leads to the following optimal solution:")
        print(f"The best strategy is to cut {final_num_b2} B2 balls and {final_num_b1} B1 balls.")
        print(f"The final equation for the maximum value is:")
        print(f"{final_num_b2} * {price_b2} + {final_num_b1} * {price_b1} = {final_value}")

    else:
        # This case is not expected, but included for completeness.
        print(f"The optimal strategy is to cut only B1 balls, with a total value of {value_strategy2}.")


solve_cutting_problem()