import math

def solve_cutting_problem():
    """
    Solves the cutting problem based on the provided formulation.
    """
    # Billet dimensions in cm
    billet_l, billet_w, billet_h = 16, 11, 4

    # Product details (price and bounding box dimension in cm)
    price_b2 = 150
    dim_b2 = 4

    price_b1 = 1
    dim_b1 = 1

    # --- Scenario 1: Using B2 and B1 pieces ---
    # Since B2s are very valuable, we prioritize placing as many as possible.
    num_b2_l = math.floor(billet_l / dim_b2)
    num_b2_w = math.floor(billet_w / dim_b2)
    num_b2_h = math.floor(billet_h / dim_b2)
    num_b2 = num_b2_l * num_b2_w * num_b2_h
    value_b2 = num_b2 * price_b2

    # The B2 pieces occupy a rectangular space of (4*4) x (2*4) x (1*4) = 16x8x4 cm.
    # The remaining space is a block of 16 x (11-8) x 4 = 16x3x4 cm.
    # As per the analysis of the constraints, this remaining space can be filled with B1 pieces.
    rem_l = billet_l
    rem_w = billet_w - (num_b2_w * dim_b2)
    rem_h = billet_h
    
    num_b1_scenario1 = math.floor(rem_l / dim_b1) * math.floor(rem_w / dim_b1) * math.floor(rem_h / dim_b1)
    value_b1_scenario1 = num_b1_scenario1 * price_b1

    total_value_scenario1 = value_b2 + value_b1_scenario1

    # --- Scenario 2: Using T1 and B1 pieces (no B2s) ---
    # The provided constraints for T1 pieces are extremely restrictive, making them
    # far less space-efficient than B1 pieces. The optimal strategy in this
    # scenario is to ignore T1s and fill the entire billet with B1 pieces.
    num_b1_scenario2 = math.floor(billet_l / dim_b1) * math.floor(billet_w / dim_b1) * math.floor(billet_h / dim_b1)
    total_value_scenario2 = num_b1_scenario2 * price_b1
    
    # --- Conclusion ---
    # Compare the total values from both scenarios to find the maximum.
    if total_value_scenario1 > total_value_scenario2:
        print(f"{num_b2} * {price_b2} + {num_b1_scenario1} * {price_b1} = {total_value_scenario1}")
    else:
        print(f"{num_b1_scenario2} * {price_b1} = {total_value_scenario2}")

solve_cutting_problem()