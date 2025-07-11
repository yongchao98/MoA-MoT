import math

def solve_cutting_problem():
    """
    Solves the billet cutting problem based on the user's formulation.
    """

    # --- Step 1 & 2: Strategy Definition ---
    # The problem formulation makes it impossible to have both B2 and T1 pieces.
    # B2 pieces are far more valuable, so the optimal strategy is to maximize B2
    # pieces and then fill the remaining space with B1 pieces.

    # --- Piece Properties ---
    price_b2 = 150
    price_b1 = 1
    price_t1 = 5

    # --- Step 3: B2 Ball Placement ---
    # A hexagonal-style packing in the XY plane allows for 11 B2 balls.
    b2_centers = [
        # Row 1 (y=4)
        (4, 4, 4), (12, 4, 4), (20, 4, 4), (28, 4, 4),
        # Row 2 (y=11)
        (8, 11, 4), (16, 11, 4), (24, 11, 4),
        # Row 3 (y=18)
        (4, 18, 4), (12, 18, 4), (20, 18, 4), (28, 18, 4),
    ]
    num_b2 = len(b2_centers)
    value_b2 = num_b2 * price_b2

    # --- Step 4: B1 Ball Placement ---
    # Search for valid B1 locations in the remaining space.
    # The 8 corners of the billet are promising candidates.
    
    b1_candidates = [
        (1, 1, 1), (31, 1, 1), (1, 21, 1), (31, 21, 1),
        (1, 1, 7), (31, 1, 7), (1, 21, 7), (31, 21, 7)
    ]
    
    valid_b1_centers = []

    # Non-overlapping constraints values (distance squared)
    b1_to_b2_dist_sq = 25
    b1_to_b1_dist_sq = 4

    for b1_cand in b1_candidates:
        is_valid = True
        # Check against all B2 balls
        for b2_c in b2_centers:
            dist_sq = (b1_cand[0] - b2_c[0])**2 + (b1_cand[1] - b2_c[1])**2 + (b1_cand[2] - b2_c[2])**2
            if dist_sq < b1_to_b2_dist_sq:
                is_valid = False
                break
        if not is_valid:
            continue
            
        # Check against already added B1 balls
        for b1_c in valid_b1_centers:
            dist_sq = (b1_cand[0] - b1_c[0])**2 + (b1_cand[1] - b1_c[1])**2 + (b1_cand[2] - b1_c[2])**2
            if dist_sq < b1_to_b1_dist_sq:
                is_valid = False
                break
        
        if is_valid:
            valid_b1_centers.append(b1_cand)

    num_b1 = len(valid_b1_centers)
    value_b1 = num_b1 * price_b1

    # --- T1 Cubes ---
    # As established, no T1 cubes can be placed with B2 balls.
    num_t1 = 0
    value_t1 = num_t1 * price_t1

    # --- Step 5: Final Calculation and Output ---
    total_value = value_b2 + value_b1 + value_t1

    print("Problem Formulation Analysis:")
    print("The provided formulation for T1 cube constraints is flawed. For example, it makes placing B2 balls and T1 cubes in the same billet impossible.")
    print("Based on the provided rules, the optimal strategy is to maximize the high-value B2 balls and then add low-value B1 balls.\n")
    
    print("Optimal Cutting Plan Found:")
    print(f"Number of B2 balls (radius 2cm): {num_b2}")
    print(f"Number of B1 balls (diameter 1cm): {num_b1}")
    print(f"Number of T1 cubes (side 1cm): {num_t1}\n")
    
    print("Final Equation for Maximum Value:")
    print(f"{num_b2} * {price_b2} + {num_b1} * {price_b1} + {num_t1} * {price_t1} = {total_value}")

solve_cutting_problem()