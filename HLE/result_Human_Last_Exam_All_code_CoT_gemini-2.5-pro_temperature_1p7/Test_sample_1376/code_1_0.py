import math

def solve_emperor_tomb_engraving():
    """
    Calculates the optimal number of squares and circles to maximize the number of engraved characters.
    """
    # Step 1: Define constants from the problem description
    sheet_w, sheet_h = 140, 110
    square_side = 10
    circle_radius = 20
    circle_diameter = circle_radius * 2

    chars_per_square = 4
    
    # Step 2: Calculate characters per circle
    # To encode 7*7*9 = 441 unique characters with 8 bagua symbols, we need 3 symbols per character (8^2=64, 8^3=512).
    # A circle holds 999 symbols, so it fits 999 // 3 characters.
    chars_per_circle = 999 // 3
    
    max_k = 0
    best_n = 0
    best_m = 0
    
    print("Analyzing the optimal cutting plan for the 140x110cm meteorite sheet...")
    print("-" * 30)

    # Helper function to calculate the number of squares in an L-shaped area
    # This is created when a rectangular block is cut from a corner.
    def calculate_n_for_l_shape(rect_w, rect_h, block_w, block_h):
        if block_w > rect_w or block_h > rect_h:
            return 0
        
        # Calculate squares in the two rectangles that form the L-shape
        squares_in_rect1 = (block_w // square_side) * ((rect_h - block_h) // square_side)
        squares_in_rect2 = ((rect_w - block_w) // square_side) * (rect_h // square_side)
        return squares_in_rect1 + squares_in_rect2

    # Step 3: Iterate through all possible numbers of circles (M), from 0 to 6.
    max_m = (sheet_w // circle_diameter) * (sheet_h // circle_diameter)
    for m in range(max_m + 1):
        n = 0
        
        # Step 4: For each M, determine the optimal packing and calculate the resulting N
        if m == 0:
            # Entire sheet is used for squares
            n = (sheet_w // square_side) * (sheet_h // square_side)
        elif m == 1:
            # 1x1 block of circles (40x40)
            n = calculate_n_for_l_shape(sheet_w, sheet_h, 1 * circle_diameter, 1 * circle_diameter)
        elif m == 2:
            # 2x1 block (80x40) vs 1x2 block (40x80)
            n1 = calculate_n_for_l_shape(sheet_w, sheet_h, 2 * circle_diameter, 1 * circle_diameter)
            n2 = calculate_n_for_l_shape(sheet_w, sheet_h, 1 * circle_diameter, 2 * circle_diameter)
            n = max(n1, n2)
        elif m == 3:
            # 3x1 block (120x40). 1x3 is not possible as 3*40 > 110.
            n = calculate_n_for_l_shape(sheet_w, sheet_h, 3 * circle_diameter, 1 * circle_diameter)
        elif m == 4:
            # 2x2 block (80x80)
            n = calculate_n_for_l_shape(sheet_w, sheet_h, 2 * circle_diameter, 2 * circle_diameter)
        elif m == 5:
            # This forms a 3x2 block (120x80) with a 40x40 hole.
            # Calculate squares in the L-shape around the 120x80 bounding box.
            n_l_shape = calculate_n_for_l_shape(sheet_w, sheet_h, 3 * circle_diameter, 2 * circle_diameter)
            # Add squares from the 40x40 empty hole inside the circle block.
            n_hole = (circle_diameter // square_side) ** 2
            n = n_l_shape + n_hole
        elif m == 6:
            # 3x2 block (120x80)
            n = calculate_n_for_l_shape(sheet_w, sheet_h, 3 * circle_diameter, 2 * circle_diameter)

        # Step 5: Calculate total characters K and check if it's a new maximum
        k = (chars_per_square * n) + (chars_per_circle * m)
        
        print(f"Case M={m}: Can fit {n} squares. Total characters K = {k}")

        if k > max_k:
            max_k = k
            best_n = n
            best_m = m

    # Step 6: Print the final optimal result
    print("-" * 30)
    print("Conclusion: The optimal strategy is to produce artifacts for the maximum character count.")
    print(f"The maximum is achieved with N = {best_n} squares and M = {best_m} circles.")
    
    k_final_term1 = chars_per_square * best_n
    k_final_term2 = chars_per_circle * best_m
    print(f"The final calculation is: {chars_per_square} * {best_n} + {chars_per_circle} * {best_m} = {k_final_term1} + {k_final_term2} = {max_k}")
    
    final_answer = f"{best_n}:{best_m}:{max_k}"
    print(f"\nThe final answer in N:M:K format is: {final_answer}")
    print(f"<<<{final_answer}>>>")

solve_emperor_tomb_engraving()