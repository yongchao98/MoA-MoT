import math

def solve_cutting_problem():
    """
    Solves the cutting problem based on the user's formulation.

    The key insight is that the T1-B2 non-overlapping constraint makes it
    impossible for T1 and B2 pieces to coexist in the same billet.
    Therefore, we only need to compare the maximum value achievable with
    (B2 + B1) against the maximum value achievable with (T1 + B1).
    """

    # --- Scenario A: B2 and B1 pieces ---
    # The value is dominated by B2 pieces. A simple grid packing allows for
    # 8 B2 pieces.
    num_b2 = 8
    price_b2 = 150
    value_scenario_a = num_b2 * price_b2 # 8 * 150 = 1200

    # --- Scenario B: T1 and B1 pieces ---
    # The value is dominated by T1 pieces. The non-overlapping constraint
    # min(|xi-xj|,...) >= 2 allows placing T1 centers on a grid with step 2.
    price_t1 = 5
    t1_center_x_range = (1, 31)
    t1_center_y_range = (1, 21)
    t1_center_z_range = (1, 7)
    step = 2
    
    # Calculate the number of possible T1s
    num_t1_x = (t1_center_x_range[1] - t1_center_x_range[0]) // step + 1
    num_t1_y = (t1_center_y_range[1] - t1_center_y_range[0]) // step + 1
    num_t1_z = (t1_center_z_range[1] - t1_center_z_range[0]) // step + 1
    num_t1 = num_t1_x * num_t1_y * num_t1_z

    # With this packing, no B1s can be added due to T1-B1 constraint.
    num_b1_in_b = 0
    price_b1 = 1
    
    value_scenario_b = num_t1 * price_t1 + num_b1_in_b * price_b1

    # --- Conclusion ---
    # Compare the scenarios and print the result for the best one.
    if value_scenario_b > value_scenario_a:
        # The best solution is from Scenario B
        final_value = value_scenario_b
        equation = f"{num_t1} * {price_t1}"
        print(f"{equation} = {final_value}")
    else:
        # The best solution is from Scenario A
        final_value = value_scenario_a
        equation = f"{num_b2} * {price_b2}"
        # Simplified: we ignore B1s as they don't change the outcome
        print(f"{equation} = {final_value}")

solve_cutting_problem()