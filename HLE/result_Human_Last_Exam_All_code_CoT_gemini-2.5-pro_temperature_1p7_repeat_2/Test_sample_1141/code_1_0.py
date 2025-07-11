import math

def solve_cutting_problem():
    """
    Calculates the maximum value based on the user's problem formulation.

    The strategy is to place the maximum number of high-value B2 balls first,
    then fill the remaining space with B1 balls, adhering to all constraints.
    As explained, strategies involving T1 cubes are not competitive due to the
    problem's specific non-overlapping constraints.
    """
    # Define prices
    price_b2 = 150
    price_b1 = 1
    price_t1 = 5

    # --- B2 Placement ---
    # Based on a grid packing, we can fit 8 B2 balls.
    b2_centers = []
    b2_x_coords = [4, 12, 20, 28]
    b2_y_coords = [4, 12]
    b2_z_coord = 4
    for x in b2_x_coords:
        for y in b2_y_coords:
            b2_centers.append((x, y, b2_z_coord))
    
    num_b2 = len(b2_centers)

    # --- B1 Placement ---
    # We create a grid of candidate centers for B1 balls. This grid automatically
    # satisfies the B1-to-B1 non-overlapping constraint:
    # (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 4
    # because any two points on the grid are separated by at least 2 units
    # along one or more axes.
    b1_candidate_centers = []
    # Valid B1 center ranges: x in [1, 31], y in [1, 21], z in [1, 7]
    for x in range(1, 32, 2):
        for y in range(1, 22, 2):
            for z in range(1, 8, 2):
                b1_candidate_centers.append((x, y, z))

    # Filter B1 candidates based on the B1-to-B2 constraint:
    # (xi-xj)^2 + (yi-yj)^2 + (zi-zj)^2 >= 25
    valid_b1_count = 0
    for b1_c in b1_candidate_centers:
        is_valid = True
        for b2_c in b2_centers:
            dist_sq = (b1_c[0] - b2_c[0])**2 + (b1_c[1] - b2_c[1])**2 + (b1_c[2] - b2_c[2])**2
            if dist_sq < 25:
                is_valid = False
                break
        if is_valid:
            valid_b1_count += 1
    
    num_b1 = valid_b1_count
    num_t1 = 0  # As explained, T1 pieces are not viable in the optimal solution.

    # --- Calculate Final Value ---
    total_value = (num_b2 * price_b2) + (num_b1 * price_b1) + (num_t1 * price_t1)

    # --- Print Final Equation ---
    print(f"{num_b2} pieces of B2 at price {price_b2} + {num_b1} pieces of B1 at price {price_b1} + {num_t1} pieces of T1 at price {price_t1} = {total_value}")
    
    # Return the answer in the requested format
    print(f"\n<<<{total_value}>>>")

solve_cutting_problem()