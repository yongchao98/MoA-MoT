import math

def solve_engraving_maximization():
    """
    Calculates the optimal number of squares (N) and circles (M) to maximize
    the total number of engraved characters (K).
    """

    # --- Step 1: Define Constants from the problem statement ---
    RECT_W, RECT_H = 140, 110
    SQUARE_SIDE = 10
    CIRCLE_DIAMETER = 40  # 20cm radius * 2

    CHARS_PER_SQUARE = 4
    SYMBOLS_PER_CIRCLE = 999
    
    # --- Step 2: Determine Characters per Circle based on problem hints ---
    # The problem's emphasis on the number 9 (7x7x9 chars, 999 symbols)
    # suggests a direct relationship. We assume 9 symbols encode 1 character.
    SYMBOLS_PER_BIO_CHAR = 9
    CHARS_PER_CIRCLE = SYMBOLS_PER_CIRCLE // SYMBOLS_PER_BIO_CHAR

    # --- Step 3: Initialize variables to store the best result ---
    best_N, best_M, max_K = 0, 0, 0

    # --- Step 4: Systematically check all possible grid layouts of circles ---
    # We check both orientations of the main rectangle (140x110 and 110x140)
    for rect_w, rect_h in [(RECT_W, RECT_H), (RECT_H, RECT_W)]:
        
        max_m_cols = rect_w // CIRCLE_DIAMETER
        max_m_rows = rect_h // CIRCLE_DIAMETER

        for m_cols in range(max_m_cols + 1):
            for m_rows in range(max_m_rows + 1):
                
                # Number of circles for this layout
                M = m_cols * m_rows
                
                # Area occupied by the grid of circles
                circles_w = m_cols * CIRCLE_DIAMETER
                circles_h = m_rows * CIRCLE_DIAMETER

                # Calculate the remaining area, which is L-shaped.
                # We can split this L-shape into two non-overlapping rectangles
                # to calculate the number of squares that can fit.
                
                # Area 1: The rectangle along the width of the main sheet
                rem_rect1_w = rect_w - circles_w
                rem_rect1_h = rect_h
                
                # Area 2: The remaining rectangle below the circle block
                rem_rect2_w = circles_w
                rem_rect2_h = rect_h - circles_h

                squares_in_rect1 = (rem_rect1_w // SQUARE_SIDE) * (rem_rect1_h // SQUARE_SIDE)
                squares_in_rect2 = (rem_rect2_w // SQUARE_SIDE) * (rem_rect2_h // SQUARE_SIDE)
                
                N = squares_in_rect1 + squares_in_rect2

                # --- Step 5: Calculate total characters (K) and check if it's the maximum so far ---
                K = (N * CHARS_PER_SQUARE) + (M * CHARS_PER_CIRCLE)
                
                if K > max_K:
                    max_K = K
                    best_N = N
                    best_M = M
    
    # --- Step 6: Print the final result in the specified N:M:K format ---
    # The output shows the optimal number of squares (N), circles (M),
    # and the corresponding maximum total characters engraved (K).
    print(f"To maximize the number of engraved characters, the workers should produce:")
    print(f"{best_N} squares and {best_M} circles.")
    print(f"This results in a maximum of {max_K} characters.")
    print("\nFinal Answer Format:")
    print(f"{best_N}:{best_M}:{max_K}")

solve_engraving_maximization()