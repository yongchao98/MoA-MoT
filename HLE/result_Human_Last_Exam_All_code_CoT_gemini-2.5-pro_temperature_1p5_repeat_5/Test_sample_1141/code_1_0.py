import math

def solve_cutting_problem():
    """
    Solves the cutting problem by comparing two mutually exclusive scenarios
    and calculating the maximum value.
    """

    # --- Problem Definition ---
    # Prices
    price_b2 = 150
    price_t1 = 5
    price_b1 = 1

    # Billet dimensions in 0.5 cm grid units
    billet_dims = {'x': 32, 'y': 22, 'z': 8}

    # Item center constraints (within billet)
    # B2 (radius 4): x in [4, 28], y in [4, 18], z = 4
    # B1 (radius 1) / T1 (side 2): x in [1, 31], y in [1, 21], z in [1, 7]

    # --- Scenario A: Packing B2 and B1 items ---

    # 1. Pack B2 items.
    # B2 balls require z_center=4. Packing is a 2D problem.
    # B2-B2 non-overlap: distance between centers >= 8 units.
    # We can place a 4x2 rectangular grid of B2 balls.
    # X-centers: 4, 12, 20, 28
    # Y-centers: 4, 12
    num_b2 = 4 * 2

    # 2. Pack B1 items in the remaining space.
    # B1-B2 non-overlap: distance^2 >= 25.
    # The B2 centers are at y=4 and y=12. A B1 center at y>=17 satisfies
    # |y_b1 - 12| >= 5, which ensures (y_b1 - y_b2)^2 >= 25.
    # This means the entire volume with centers in [1,31]x[17,21]x[1,7] is available for B1s.
    # B1-B1 non-overlap: distance^2 >= 4. A grid with step 2 works.
    x_pos_b1 = math.floor((31 - 1) / 2) + 1  # 16 positions
    y_pos_b1 = math.floor((21 - 17) / 2) + 1 # 3 positions (17, 19, 21)
    z_pos_b1 = math.floor((7 - 1) / 2) + 1   # 4 positions (1, 3, 5, 7)
    num_b1_scenario_a = x_pos_b1 * y_pos_b1 * z_pos_b1

    value_scenario_a = num_b2 * price_b2 + num_b1_scenario_a * price_b1

    # --- Scenario B: Packing T1 and B1 items ---
    
    # Analysis showed that due to the very strict T1-B1 constraint
    # (min(|dx|,|dy|,|dz|) >= 2), adding a T1 (value 5) eliminates
    # far more B1s (value 1) than it's worth.
    # For example, one T1 at (1,1,1) eliminates over 250 B1 grid points.
    # The optimal strategy for this scenario is therefore to place zero T1s
    # and fill the entire billet with B1s.
    num_t1 = 0
    x_pos_b1_full = math.floor((31 - 1) / 2) + 1 # 16 positions
    y_pos_b1_full = math.floor((21 - 1) / 2) + 1 # 11 positions
    z_pos_b1_full = math.floor((7 - 1) / 2) + 1  # 4 positions
    num_b1_scenario_b = x_pos_b1_full * y_pos_b1_full * z_pos_b1_full

    value_scenario_b = num_t1 * price_t1 + num_b1_scenario_b * price_b1
    
    # --- Final Result ---
    # The problem asks for the equation for the highest valid solution.
    # We compare the two scenarios.
    
    if value_scenario_a > value_scenario_b:
        num_b1_final = num_b1_scenario_a
        final_value = value_scenario_a
        # Print the final equation with each number and operator separated by spaces.
        print(num_b2, "*", price_b2, "+", num_b1_final, "*", price_b1, "=", final_value)

    else:
        # This case is not expected to be chosen, but included for completeness.
        num_b1_final = num_b1_scenario_b
        final_value = value_scenario_b
        print(num_t1, "*", price_t1, "+", num_b1_final, "*", price_b1, "=", final_value)


solve_cutting_problem()
>>> 8 * 150 + 192 * 1 = 1392