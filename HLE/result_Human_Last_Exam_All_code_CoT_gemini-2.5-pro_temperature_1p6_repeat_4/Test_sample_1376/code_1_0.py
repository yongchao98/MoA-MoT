import math

def solve_engraving_optimization():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # --- Define problem constants ---
    sheet_w = 140
    sheet_h = 110

    square_size = 10
    # A circle with a 20cm radius fits in a 40x40cm square bounding box.
    circle_bbox_size = 40

    square_chars = 4
    circle_chars = 999

    print("--- Problem Analysis ---")
    print("The goal is to maximize the total characters K from N squares and M circles.")
    print(f"The optimization formula is: K = ({square_chars} * N) + ({circle_chars} * M)")
    print("The strategy is to prioritize cutting circles, as they are more 'valuable' per area.")
    print("\n--- Step 1: Maximize the Number of Circles (M) ---")

    # Calculate how many 40x40cm circle placeholders fit on the 140x110cm sheet.
    # Note: floor division (//) is used as we can only cut whole items.
    circles_along_w = sheet_w // circle_bbox_size
    circles_along_h = sheet_h // circle_bbox_size
    M = circles_along_w * circles_along_h

    print(f"Sheet dimensions: {sheet_w}cm x {sheet_h}cm")
    print(f"Circle material required: {circle_bbox_size}cm x {circle_bbox_size}cm")
    print(f"Circles that fit along the 140cm side: {circles_along_w}")
    print(f"Circles that fit along the 110cm side: {circles_along_h}")
    print(f"Maximum number of circles (M) = {circles_along_w} * {circles_along_h} = {M}")

    # --- Step 2: Calculate Remaining Area for Squares (N) ---
    # The M circles will be cut from a solid block of material.
    circles_block_w = circles_along_w * circle_bbox_size
    circles_block_h = circles_along_h * circle_bbox_size

    # This leaves two rectangular strips of material.
    # Strip 1: Along the full width of the sheet.
    rem_strip1_w = sheet_w
    rem_strip1_h = sheet_h - circles_block_h
    # Strip 2: Beside the block of circles.
    rem_strip2_w = sheet_w - circles_block_w
    rem_strip2_h = circles_block_h

    print("\n--- Step 2: Maximize the Number of Squares (N) in Remaining Area ---")
    print(f"The {M} circles occupy a {circles_block_w}x{circles_block_h}cm area.")
    print(f"This leaves two strips of material: {rem_strip1_w}x{rem_strip1_h}cm and {rem_strip2_w}x{rem_strip2_h}cm.")
    
    squares_in_strip1 = (rem_strip1_w // square_size) * (rem_strip1_h // square_size)
    squares_in_strip2 = (rem_strip2_w // square_size) * (rem_strip2_h // square_size)
    N = squares_in_strip1 + squares_in_strip2
    
    print(f"Squares from first strip: ({rem_strip1_w}//{square_size}) * ({rem_strip1_h}//{square_size}) = {squares_in_strip1}")
    print(f"Squares from second strip: ({rem_strip2_w}//{square_size}) * ({rem_strip2_h}//{square_size}) = {squares_in_strip2}")
    print(f"Total number of squares (N) = {squares_in_strip1} + {squares_in_strip2} = {N}")
    
    # --- Step 3: Calculate the Maximum Characters (K) ---
    chars_from_squares = square_chars * N
    chars_from_circles = circle_chars * M
    K = chars_from_squares + chars_from_circles
    
    print("\n--- Final Result ---")
    print(f"Optimal number of squares (N): {N}")
    print(f"Optimal number of circles (M): {M}")
    print("\nThe maximum number of characters K is found with the equation:")
    print(f"K = ({square_chars} * {N}) + ({circle_chars} * {M}) = {chars_from_squares} + {chars_from_circles} = {K}")
    
    print(f"\nFinal answer in N:M:K format is: {N}:{M}:{K}")

solve_engraving_optimization()