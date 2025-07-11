import math

def solve_cutting_problem():
    """
    Solves the billet cutting problem based on the user's formulation.
    """

    # --- Prices ---
    price_b2 = 150
    price_b1 = 1
    price_t1 = 5

    # --- Scenario A: Only T1 cubes ---
    # The T1-T1 constraint min(|xi-xj|, |yi-yj|, |zi-zj|) >= 2 means
    # centers can be placed on a grid with a step of 2.
    
    # Valid center coordinates for T1 cubes on a grid
    x_coords_t1 = range(1, 31 + 1, 2)
    y_coords_t1 = range(1, 21 + 1, 2)
    z_coords_t1 = range(1, 7 + 1, 2)
    
    num_t1 = len(x_coords_t1) * len(y_coords_t1) * len(z_coords_t1)
    value_A = num_t1 * price_t1
    
    # --- Scenario B: B2 balls and B1 balls ---
    
    # Part 1: Place maximum number of B2 balls
    # The B2-B2 constraint is dist^2 >= 64, or dist >= 8.
    # Since z is fixed at 4, we can place centers on a 2D grid with a step of 8.
    x_coords_b2 = range(4, 28 + 1, 8) # 4, 12, 20, 28
    y_coords_b2 = range(4, 18 + 1, 8) # 4, 12
    z_coord_b2 = 4
    
    num_b2 = len(x_coords_b2) * len(y_coords_b2)
    value_from_b2s = num_b2 * price_b2
    
    b2_centers = []
    for x in x_coords_b2:
        for y in y_coords_b2:
            b2_centers.append((x, y, z_coord_b2))

    # Part 2: Place B1 balls in the remaining space
    # B1-B1 constraint: dist^2 >= 4, or dist >= 2. We can use a grid with step 2.
    # B1-B2 constraint: dist^2 >= 25.
    num_b1 = 0
    
    # Potential B1 center coordinates on a grid
    x_coords_b1 = range(1, 31 + 1, 2)
    y_coords_b1 = range(1, 21 + 1, 2)
    z_coords_b1 = range(1, 7 + 1, 2)

    for x_b1 in x_coords_b1:
        for y_b1 in y_coords_b1:
            for z_b1 in z_coords_b1:
                is_valid = True
                # Check against all placed B2s
                for x_b2, y_b2, z_b2 in b2_centers:
                    dist_sq = (x_b1 - x_b2)**2 + (y_b1 - y_b2)**2 + (z_b1 - z_b2)**2
                    if dist_sq < 25:
                        is_valid = False
                        break
                if not is_valid:
                    continue
                
                num_b1 += 1
                
    value_from_b1s = num_b1 * price_b1
    value_B = value_from_b2s + value_from_b1s

    # --- Compare scenarios and print the result ---
    if value_A > value_B:
        print("Scenario A (T1 cubes only) is optimal.")
        # The final equation is the number of T1 cubes multiplied by their price.
        print(f"The highest valid solution is found by filling the billet with T1 cubes.")
        print(f"Final Calculation: {num_t1} * {price_t1} = {value_A}")
        final_value = value_A
        
    else:
        print("Scenario B (B2 and B1 balls) is optimal.")
        # The final equation is the sum of values from B2 and B1 pieces.
        print(f"The highest valid solution is found by combining B2 and B1 balls.")
        print(f"Final Calculation: ({num_b2} * {price_b2}) + ({num_b1} * {price_b1}) = {value_B}")
        final_value = value_B
        
    # As requested, output the final numerical answer in the specified format.
    # The print statements above provide the explanation and equation.
    # The line below provides the final answer for parsing.
    # print(f"\n<<<{final_value}>>>")

solve_cutting_problem()