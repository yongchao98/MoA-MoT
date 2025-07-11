import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize the
    number of engraved characters on a given sheet of material.
    """
    # --- Problem Parameters ---
    MATERIAL_W = 140  # cm
    MATERIAL_H = 110  # cm
    SQUARE_SIDE = 10  # cm
    CIRCLE_DIAMETER = 40  # cm (from 20cm radius)
    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    print("--- Step 1: Maximizing the number of circles (M) ---")
    # The character density of circles is much higher, so we prioritize them.
    # We need to fit 40x40cm bounding boxes into the 140x110cm material.
    
    # Calculate how many circles fit along each dimension
    fit_along_w = MATERIAL_W // CIRCLE_DIAMETER
    fit_along_h = MATERIAL_H // CIRCLE_DIAMETER
    
    M = fit_along_w * fit_along_h
    
    print(f"The material is {MATERIAL_W}x{MATERIAL_H} cm. A circle requires a {CIRCLE_DIAMETER}x{CIRCLE_DIAMETER} cm area.")
    print(f"Fit along width: floor({MATERIAL_W} / {CIRCLE_DIAMETER}) = {fit_along_w}")
    print(f"Fit along height: floor({MATERIAL_H} / {CIRCLE_DIAMETER}) = {fit_along_h}")
    print(f"Total circles (M) = {fit_along_w} * {fit_along_h} = {M}")
    print("-" * 20)

    print("--- Step 2: Calculating squares (N) in the remaining area ---")
    # The 6 circles are arranged in a 3x2 grid, occupying an area.
    circles_area_w = fit_along_w * CIRCLE_DIAMETER
    circles_area_h = fit_along_h * CIRCLE_DIAMETER
    print(f"The {M} circles occupy a {circles_area_w}x{circles_area_h} cm area.")

    # The remaining area is an L-shape, which can be split into two rectangles.
    # Rectangle A: (MATERIAL_W - circles_area_w) x MATERIAL_H
    rem_rect_A_w = MATERIAL_W - circles_area_w
    rem_rect_A_h = MATERIAL_H
    
    # Rectangle B: circles_area_w x (MATERIAL_H - circles_area_h)
    rem_rect_B_w = circles_area_w
    rem_rect_B_h = MATERIAL_H - circles_area_h

    print(f"The remaining space can be split into two rectangles:")
    print(f" - Rectangle A: {rem_rect_A_w}x{rem_rect_A_h} cm")
    print(f" - Rectangle B: {rem_rect_B_w}x{rem_rect_B_h} cm")

    # Calculate how many 10x10 squares fit in each remaining rectangle.
    squares_in_A = (rem_rect_A_w // SQUARE_SIDE) * (rem_rect_A_h // SQUARE_SIDE)
    squares_in_B = (rem_rect_B_w // SQUARE_SIDE) * (rem_rect_B_h // SQUARE_SIDE)
    
    N = squares_in_A + squares_in_B
    
    print(f"Squares in Rectangle A = ({rem_rect_A_w}/{SQUARE_SIDE}) * ({rem_rect_A_h}/{SQUARE_SIDE}) = {squares_in_A}")
    print(f"Squares in Rectangle B = ({rem_rect_B_w}/{SQUARE_SIDE}) * ({rem_rect_B_h}/{SQUARE_SIDE}) = {squares_in_B}")
    print(f"Total squares (N) = {squares_in_A} + {squares_in_B} = {N}")
    print("-" * 20)

    print("--- Step 3: Calculating the maximum total characters (K) ---")
    # Calculate K using the formula K = 4*N + 999*M
    K = (CHARS_PER_SQUARE * N) + (CHARS_PER_CIRCLE * M)
    
    print(f"The formula for total characters is: K = (chars_per_square * N) + (chars_per_circle * M)")
    print(f"K = ({CHARS_PER_SQUARE} * {N}) + ({CHARS_PER_CIRCLE} * {M})")
    print(f"K = {CHARS_PER_SQUARE * N} + {CHARS_PER_CIRCLE * M} = {K}")
    print("-" * 20)

    print("--- Final Answer ---")
    print("The optimal production is N:M:K")
    print(f"{N}:{M}:{K}")

solve_engraving_problem()