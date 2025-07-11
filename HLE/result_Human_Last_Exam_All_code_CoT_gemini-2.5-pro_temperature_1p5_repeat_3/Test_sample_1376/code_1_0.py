import math

def solve_engraving_problem():
    """
    Calculates the maximum number of characters that can be engraved
    by finding the optimal number of squares (N) and circles (M).
    """
    # Define the problem's parameters
    RECT_W = 140  # cm
    RECT_H = 110  # cm
    SQUARE_SIDE = 10  # cm
    CIRCLE_DIAMETER = 40  # cm
    CIRCLE_RADIUS = 20 # cm
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    # The core of the problem is to maximize K = N * 4 + M * 999.
    # Given the high value of characters per circle, we must maximize M.
    # We investigate hexagonal packing, which is denser than simple grid packing.

    # --- Analysis of Hexagonal Packing ---
    # We determine if we can fit more circles than the simple grid packing (which allows 3x2=6 circles).
    # Let's arrange circles with rows parallel to the 140cm side.
    # A staggered arrangement allows for 3 rows:
    # Row 1 centers (y=20): x=20, 60, 100 (3 circles)
    # Row 2 centers (y~=54.64): x=40, 80, 120 (3 circles)
    # Row 3 centers (y~=89.28): x=20, 60, 100 (3 circles)
    m_hex = 9

    # We calculate the bounding box of this 9-circle arrangement.
    # Max x-center is 120, max y-center is ~89.28
    bbox_hex_w = 120 + CIRCLE_RADIUS  # 120 + 20 = 140 cm
    # The height is determined by the y-position of the center of the last row plus the radius.
    # y_pitch = r * sqrt(3) ~= 34.64
    # y_center_last_row = r + 2 * y_pitch ~= 20 + 69.28 = 89.28
    bbox_hex_h = (CIRCLE_RADIUS + 2 * (CIRCLE_RADIUS * math.sqrt(3))) + CIRCLE_RADIUS
    bbox_hex_h = 89.28427 + 20 # Approximation based on calculation
    
    # This arrangement fits perfectly in a 140cm x 109.28cm box.
    # The arrangement fits inside the 140x110cm material.
    
    # The remaining space is a thin strip of 140cm x (110 - 109.28)cm.
    # This strip is too narrow (0.72cm) to fit any 10cm squares.
    n_hex = 0

    # The hexagonal packing yields M=9 circles and N=0 squares.
    # This is superior to grid packing (M=6), as sacrificing even one circle (999 chars)
    # for squares (max gain of 16*4=64 chars) results in a net loss.
    # Therefore, this is the optimal solution.

    final_N = n_hex
    final_M = m_hex
    final_K = final_N * CHARS_PER_SQUARE + final_M * CHARS_PER_CIRCLE

    print("To maximize the number of characters, the densest packing for circles should be used.")
    print("A hexagonal packing arrangement allows for 9 circles to be cut from the 140x110cm material.")
    print("This arrangement leaves no significant space to cut any 10x10cm squares.")
    print("\n--- Optimal Solution ---")
    print(f"Number of squares (N): {final_N}")
    print(f"Number of circles (M): {final_M}")
    print(f"Maximum characters (K): {final_K}")

    print("\n--- Final Calculation ---")
    print(f"The equation for the maximum number of characters is:")
    print(f"{final_K} = {final_N} * {CHARS_PER_SQUARE} + {final_M} * {CHARS_PER_CIRCLE}")
    
    print(f"\nThe final answer in N:M:K format is: {final_N}:{final_M}:{final_K}")

solve_engraving_problem()
<<<0:9:8991>>>