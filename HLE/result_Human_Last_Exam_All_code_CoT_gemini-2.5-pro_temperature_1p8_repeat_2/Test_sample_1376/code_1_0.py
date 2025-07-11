import math

def solve_engraving_optimization():
    """
    Solves the problem of maximizing the number of engraved Chinese characters
    by finding the optimal number of squares (N) and circles (M).
    """

    # --- Problem Definition ---
    # Material and artifact dimensions
    board_w, board_h = 140, 110
    square_dim = 10
    circle_bbox_dim = 40  # A 20cm radius circle fits in a 40x40cm square

    # Value of each artifact in terms of characters
    chars_per_square = 4
    # A circle holds 999 bagua symbols. We need log8(441) ~= 2.92 symbols per character.
    # So we use 3 symbols per character.
    # Characters per circle = 999 symbols / 3 symbols/char = 333.
    chars_per_circle = 333

    print("--- Analyzing the Optimization Problem ---")
    print(f"Goal: Maximize K = {chars_per_square}*N + {chars_per_circle}*M")
    print(f"Constraints: Packing {square_dim}x{square_dim}cm squares and {circle_bbox_dim}x{circle_bbox_dim}cm circles on a {board_w}x{board_h}cm sheet.\n")

    # --- Packing Strategy ---
    # We prioritize circles. Let's find the max possible number of circles in a grid.
    max_circles_x = board_w // circle_bbox_dim
    max_circles_y = board_h // circle_bbox_dim
    max_m = max_circles_x * max_circles_y

    # This grid of circles occupies a certain area.
    grid_w = max_circles_x * circle_bbox_dim
    grid_h = max_circles_y * circle_bbox_dim
    
    # We calculate the number of squares that fit in the remaining L-shaped area when the grid is full (M=6).
    # The L-shape is split into two rectangles:
    # 1. (board_w - grid_w) x board_h
    # 2. grid_w x (board_h - grid_h)
    rem_area1_w = board_w - grid_w
    rem_area1_h = board_h
    rem_area2_w = grid_w
    rem_area2_h = board_h - grid_h

    squares_in_rem_area1 = (rem_area1_w // square_dim) * (rem_area1_h // square_dim)
    squares_in_rem_area2 = (rem_area2_w // square_dim) * (rem_area2_h // square_dim)
    
    # This is the base number of squares if all 6 circle slots are taken by circles.
    n_base_for_max_m = squares_in_rem_area1 + squares_in_rem_area2
    
    # A single circle slot, if not used for a circle, can fit squares.
    squares_per_circle_slot = (circle_bbox_dim // square_dim) ** 2

    # --- Calculation ---
    max_k = 0
    best_n, best_m = 0, 0
    
    print("--- Iterating Through Possible Scenarios (M=0 to 6) ---")

    for m_count in range(max_m + 1):
        # For a given number of circles M, the remaining slots in the grid can be used for squares.
        empty_slots = max_m - m_count
        n_count = n_base_for_max_m + (empty_slots * squares_per_circle_slot)
        
        # Calculate the total characters K
        k_total = (chars_per_square * n_count) + (chars_per_circle * m_count)
        
        # Keep track of the maximum K found so far
        if k_total > max_k:
            max_k = k_total
            best_n = n_count
            best_m = m_count

    # --- Final Result ---
    print("\n--- Optimal Solution Found ---")
    print(f"The analysis shows that the value of K = {chars_per_square}*N + {chars_per_circle}*M increases with each additional circle.")
    print("Therefore, the maximum K is achieved when we produce the maximum possible number of circles.")
    print(f"\nThe optimal configuration is N={best_n} squares and M={best_m} circles.")
    print("The final equation for the maximal number of characters is:")
    print(f"K = {chars_per_square} * {best_n} + {chars_per_circle} * {best_m} = {max_k}")
    
    # Final answer in the required format
    # print(f"\n<<<{best_n}:{best_m}:{max_k}>>>")

solve_engraving_optimization()