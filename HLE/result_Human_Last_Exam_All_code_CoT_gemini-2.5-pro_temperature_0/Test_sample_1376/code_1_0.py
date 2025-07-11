import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize
    the total number of engraved characters.
    """
    # --- Problem Parameters ---
    RECT_W = 140  # cm
    RECT_H = 110  # cm
    SQUARE_SIDE = 10  # cm
    CIRCLE_RADIUS = 20  # cm
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    # --- Calculations ---

    # Step 1: Determine the geometry of the items to be packed
    circle_diameter = CIRCLE_RADIUS * 2

    print("To maximize the number of engraved characters, we must find the optimal number of circles and squares.")
    print(f"The material is a {RECT_W}x{RECT_H}cm rectangle.")
    print(f"Circles require a {int(circle_diameter)}x{int(circle_diameter)}cm area, and squares are {int(SQUARE_SIDE)}x{int(SQUARE_SIDE)}cm.")

    # Step 2: Maximize the number of circles (M)
    # We assume a simple grid packing, as circles provide the most value.
    circles_along_w = math.floor(RECT_W / circle_diameter)
    circles_along_h = math.floor(RECT_H / circle_diameter)
    M = circles_along_w * circles_along_h

    print("\nStep 1: Maximize the number of circles (M).")
    print(f"We can fit {circles_along_w} circles along the {RECT_W}cm side and {circles_along_h} along the {RECT_H}cm side.")
    print(f"This gives a total of M = {circles_along_w} * {circles_along_h} = {M} circles.")
    
    # The area occupied by the grid of circles
    occupied_w = circles_along_w * circle_diameter
    occupied_h = circles_along_h * circle_diameter
    print(f"These circles occupy a {int(occupied_w)}x{int(occupied_h)}cm area.")

    # Step 3: Maximize the number of squares (N) in the remaining area.
    # The remaining area is an L-shape, which can be split into two rectangles.
    # Rectangle 1:
    rem_strip1_w = RECT_W - occupied_w
    rem_strip1_h = RECT_H
    squares_in_strip1 = math.floor(rem_strip1_w / SQUARE_SIDE) * math.floor(rem_strip1_h / SQUARE_SIDE)

    # Rectangle 2:
    rem_strip2_w = occupied_w
    rem_strip2_h = RECT_H - occupied_h
    squares_in_strip2 = math.floor(rem_strip2_w / SQUARE_SIDE) * math.floor(rem_strip2_h / SQUARE_SIDE)

    N = squares_in_strip1 + squares_in_strip2
    
    print("\nStep 2: Maximize the number of squares (N) in the remaining area.")
    print("The remaining L-shaped area can be divided into two rectangles:")
    print(f" - A {int(rem_strip1_w)}x{int(rem_strip1_h)}cm rectangle, which fits {squares_in_strip1} squares.")
    print(f" - A {int(rem_strip2_w)}x{int(rem_strip2_h)}cm rectangle, which fits {squares_in_strip2} squares.")
    print(f"This gives a total of N = {squares_in_strip1} + {squares_in_strip2} = {N} squares.")

    # Step 4: Calculate the total number of characters (K)
    K = N * CHARS_PER_SQUARE + M * CHARS_PER_CIRCLE
    
    print("\nStep 3: Calculate the total number of characters (K).")
    print("Each square has 4 characters, and each circle has 999 characters.")
    print(f"The final equation is: K = (N * {CHARS_PER_SQUARE}) + (M * {CHARS_PER_CIRCLE})")
    print(f"K = ({N} * {CHARS_PER_SQUARE}) + ({M} * {CHARS_PER_CIRCLE})")
    print(f"K = {N * CHARS_PER_SQUARE} + {M * CHARS_PER_CIRCLE} = {K}")

    print("\n---")
    print("Final Answer:")
    print(f"The optimal production is N={N} squares and M={M} circles, for a total of K={K} characters.")
    print(f"The result in the format N:M:K is: {N}:{M}:{K}")

solve_engraving_problem()