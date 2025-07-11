import math

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize engraved characters.
    """
    # 1. Define constants from the problem
    MATERIAL_W, MATERIAL_H = 140, 110
    SQUARE_SIZE = 10
    CIRCLE_RADIUS = 20
    CIRCLE_BBOX_SIZE = CIRCLE_RADIUS * 2

    CHARS_PER_SQUARE = 4
    CHARS_PER_CIRCLE = 999

    # 2. Strategy: Maximize circles (M) first
    # Calculate max circles that can fit in a grid
    circles_along_w = math.floor(MATERIAL_W / CIRCLE_BBOX_SIZE)
    circles_along_h = math.floor(MATERIAL_H / CIRCLE_BBOX_SIZE)
    M = circles_along_w * circles_along_h
    
    # 3. Calculate remaining area after cutting the block of circles
    # The block of circles occupies (circles_along_w * 40) x (circles_along_h * 40)
    circles_block_w = circles_along_w * CIRCLE_BBOX_SIZE
    circles_block_h = circles_along_h * CIRCLE_BBOX_SIZE

    # The remaining area is split into two non-overlapping rectangles
    # by placing the circle block at a corner of the material.
    rem_rect1_w = MATERIAL_W - circles_block_w
    rem_rect1_h = MATERIAL_H

    rem_rect2_w = circles_block_w
    rem_rect2_h = MATERIAL_H - circles_block_h

    # 4. Calculate max squares (N) from the remaining area
    squares_from_rect1 = math.floor(rem_rect1_w / SQUARE_SIZE) * math.floor(rem_rect1_h / SQUARE_SIZE)
    squares_from_rect2 = math.floor(rem_rect2_w / SQUARE_SIZE) * math.floor(rem_rect2_h / SQUARE_SIZE)
    
    N = squares_from_rect1 + squares_from_rect2

    # 5. Calculate total characters (K)
    chars_from_squares = N * CHARS_PER_SQUARE
    chars_from_circles = M * CHARS_PER_CIRCLE
    K = chars_from_squares + chars_from_circles

    # 6. Print the results step-by-step
    print(f"Optimal number of circles (M): {M}")
    print(f"Optimal number of squares (N): {N}")
    print("-" * 30)
    print("Calculating the maximum number of characters (K):")
    
    # Per the instructions, showing the final equation with its numbers
    print(f"Equation: K = (N * {CHARS_PER_SQUARE}) + (M * {CHARS_PER_CIRCLE})")
    print(f"Calculation: K = ({N} * {CHARS_PER_SQUARE}) + ({M} * {CHARS_PER_CIRCLE})")
    print(f"Result: K = {chars_from_squares} + {chars_from_circles}")
    print(f"Total K = {K}")
    print("-" * 30)
    
    # Final answer in the specified format
    print("The final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")


solve_engraving_problem()
<<<58:6:6226>>>