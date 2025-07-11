import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize the
    total number of engraved characters on a meteorite plate.
    """
    # Define constants based on the problem description
    PLATE_W = 140  # cm
    PLATE_H = 110  # cm
    SQUARE_SIDE = 10  # cm
    # A circle of 20cm radius fits in a 40x40cm square bounding box
    CIRCLE_BB_SIDE = 40  # cm
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    def calculate_optimal_layout(plate_w, plate_h):
        """
        Calculates N, M, and K for a given plate orientation (width, height).
        The strategy is to maximize circles first, then fill the rest with squares.
        """
        # --- Step 1: Maximize the number of circles (M) ---
        # Pack circles in a grid to find the maximum number that can fit.
        cols_for_circles = plate_w // CIRCLE_BB_SIDE
        rows_for_circles = plate_h // CIRCLE_BB_SIDE
        num_circles = cols_for_circles * rows_for_circles

        # --- Step 2: Calculate space remaining for squares (N) ---
        # The area used by the grid of circles
        width_used_by_circles = cols_for_circles * CIRCLE_BB_SIDE
        height_used_by_circles = rows_for_circles * CIRCLE_BB_SIDE

        # The remaining area forms an L-shape, which can be seen as two rectangles.
        # Rectangle A (the strip along the remaining width)
        rem_rect_A_w = plate_w - width_used_by_circles
        rem_rect_A_h = plate_h
        squares_in_A = (rem_rect_A_w // SQUARE_SIDE) * (rem_rect_A_h // SQUARE_SIDE)

        # Rectangle B (the strip along the remaining height, under the circles)
        rem_rect_B_w = width_used_by_circles
        rem_rect_B_h = plate_h - height_used_by_circles
        squares_in_B = (rem_rect_B_w // SQUARE_SIDE) * (rem_rect_B_h // SQUARE_SIDE)

        num_squares = squares_in_A + squares_in_B

        # --- Step 3: Calculate the total number of characters (K) ---
        total_chars = (num_squares * CHARS_PER_SQUARE) + (num_circles * CHARS_PER_CIRCLE)

        return num_squares, num_circles, total_chars

    # Calculate for the first orientation: 140x110
    n1, m1, k1 = calculate_optimal_layout(PLATE_W, PLATE_H)

    # Calculate for the second orientation: 110x140
    n2, m2, k2 = calculate_optimal_layout(PLATE_H, PLATE_W)

    # Determine the best result between the two orientations
    if k1 >= k2:
        N, M, K = n1, m1, k1
    else:
        N, M, K = n2, m2, k2

    # --- Step 4: Print the final results ---
    print("To maximize the number of engraved characters, the workers should produce:")
    print(f"- {M} circles for the bio.")
    print(f"- {N} squares for the official name.")
    print("\nThe calculation for the maximum number of characters (K) is:")
    print(f"K = (Characters per Square * N) + (Characters per Circle * M)")
    print(f"K = ({CHARS_PER_SQUARE} * {N}) + ({CHARS_PER_CIRCLE} * {M})")
    print(f"K = {CHARS_PER_SQUARE * N} + {CHARS_PER_CIRCLE * M}")
    print(f"K = {K}")
    
    print("\nThe final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")


solve_engraving_problem()