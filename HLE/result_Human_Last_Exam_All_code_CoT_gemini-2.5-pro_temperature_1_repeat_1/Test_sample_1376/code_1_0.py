import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize the
    number of engraved characters, and prints the detailed solution.
    """
    # --- Problem Parameters ---
    # Dimensions of the meteorite material rectangle
    rect_w = 140
    rect_h = 110

    # Dimensions of the items to be cut
    square_s = 10
    # A circle with 20cm radius requires a 40x40cm bounding box
    circle_d = 40

    # Engraving capacity for each item
    chars_per_square = 4
    chars_per_circle = 999

    # --- Optimization ---
    # Strategy: Maximize the number of circles first due to their high character value.
    # We use a simple grid packing for the circle bounding boxes.

    # Calculate how many circles can fit in a grid layout
    # We check both orientations of the material, but 140/40 * 110/40 is the same as 110/40 * 140/40
    m_cols = rect_w // circle_d
    m_rows = rect_h // circle_d
    M = m_cols * m_rows

    # Calculate the area occupied by the grid of circles
    circles_used_w = m_cols * circle_d
    circles_used_h = m_rows * circle_d

    # The remaining space is composed of two rectangles when the circle grid is
    # placed in a corner.
    # Area 1: The vertical strip remaining
    rem_area1_w = rect_w - circles_used_w
    rem_area1_h = rect_h
    n1 = (rem_area1_w // square_s) * (rem_area1_h // square_s)

    # Area 2: The horizontal strip remaining
    rem_area2_w = circles_used_w
    rem_area2_h = rect_h - circles_used_h
    n2 = (rem_area2_w // square_s) * (rem_area2_h // square_s)

    # Total number of squares
    N = n1 + n2

    # --- Calculate the final result ---
    K = (chars_per_square * N) + (chars_per_circle * M)

    # --- Print the step-by-step solution ---
    print(f"To maximize the total characters, we must first maximize the number of circles.")
    print(f"The material is {rect_w}x{rect_h}cm. Circles need a {circle_d}x{circle_d}cm space.")
    print(f"Maximum number of circles (M) in a grid: floor({rect_w}/{circle_d}) * floor({rect_h}/{circle_d}) = {m_cols} * {m_rows} = {M}")
    print(f"These {M} circles occupy a {circles_used_w}x{circles_used_h}cm area, leaving space for squares.")

    print(f"\nCalculating squares (N) in the remaining space:")
    print(f"1. A {rem_area1_w}x{rem_area1_h}cm area, fitting {n1} squares.")
    print(f"2. A {rem_area2_w}x{rem_area2_h}cm area, fitting {n2} squares.")
    print(f"Total squares (N) = {n1} + {n2} = {N}")

    print(f"\nOptimal configuration: N = {N} squares, M = {M} circles.")

    print("\nCalculating maximum characters (K):")
    print(f"K = (Characters per Square * N) + (Characters per Circle * M)")
    # As requested, showing each number in the final equation
    print(f"K = ({chars_per_square} * {N}) + ({chars_per_circle} * {M})")
    print(f"K = {chars_per_square * N} + {chars_per_circle * M}")
    print(f"K = {K}")

    final_answer_string = f"{N}:{M}:{K}"
    print(f"\nThe final answer in N:M:K format is: {final_answer_string}")
    print(f"\n<<<{final_answer_string}>>>")


solve_engraving_problem()