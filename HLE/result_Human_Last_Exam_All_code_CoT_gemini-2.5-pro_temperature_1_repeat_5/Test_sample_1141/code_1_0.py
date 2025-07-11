import math

def solve_cutting_problem():
    """
    This function calculates the maximum value from cutting a steel billet
    based on the user's provided (though flawed) problem formulation.

    The strategy is to prioritize the high-value B2 balls and then fill the
    remaining space with B1 balls. T1 cubes are not viable due to the
    extreme non-overlapping constraints.
    """

    # A staggered packing allows for 11 B2 balls (radius 2cm = 4 units).
    # Their centers are fixed at z=4.
    # Non-overlap for B2-B2: (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 64
    b2_centers = [
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        (8, 11, 4), (16, 11, 4), (24, 11, 4),
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4)
    ]
    num_b2 = len(b2_centers)
    price_b2 = 150
    value_b2 = num_b2 * price_b2

    # Greedily pack B1 balls (radius 0.5cm = 1 unit) into the remaining space.
    # B1 center ranges: x in [1,31], y in [1,21], z in [1,7]
    # Constraints:
    # vs B2: (x1-xb)^2 + (y1-yb)^2 + (z1-zb)^2 >= 25
    # vs B1: (x1-x2)^2 + (y1-y2)^2 + (z1-z2)^2 >= 4
    placed_b1s = []
    
    # Iterate through all possible B1 center locations in a defined order.
    for z in range(1, 8):
        for y in range(1, 22):
            for x in range(1, 32):
                candidate_b1 = (x, y, z)
                is_valid = True

                # Check against all B2 balls
                for b2 in b2_centers:
                    dist_sq = (candidate_b1[0] - b2[0])**2 + \
                              (candidate_b1[1] - b2[1])**2 + \
                              (candidate_b1[2] - b2[2])**2
                    if dist_sq < 25:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # Check against already placed B1 balls
                for b1 in placed_b1s:
                    dist_sq = (candidate_b1[0] - b1[0])**2 + \
                              (candidate_b1[1] - b1[1])**2 + \
                              (candidate_b1[2] - b1[2])**2
                    if dist_sq < 4:
                        is_valid = False
                        break
                if not is_valid:
                    continue

                # If all checks pass, place the B1 ball
                placed_b1s.append(candidate_b1)

    num_b1 = len(placed_b1s)
    price_b1 = 1
    value_b1 = num_b1 * price_b1

    # T1 cubes cannot be placed with B2 balls, so their count and value are zero.
    num_t1 = 0
    price_t1 = 5
    value_t1 = 0

    total_value = value_b2 + value_b1 + value_t1

    print("Based on a strict interpretation of your problem formulation, the highest value is achieved by placing:")
    print(f"- {num_b2} B2 balls (2cm radius)")
    print(f"- {num_b1} B1 balls (1cm diameter)")
    print(f"- {num_t1} T1 cubes (1cm side)")
    print("\nThe calculation for the total value is:")
    print(f"{num_b2} * {price_b2} + {num_b1} * {price_b1} + {num_t1} * {price_t1} = {total_value}")
    
    # The final answer is wrapped for the system.
    print(f"\n<<<{total_value}>>>")

solve_cutting_problem()