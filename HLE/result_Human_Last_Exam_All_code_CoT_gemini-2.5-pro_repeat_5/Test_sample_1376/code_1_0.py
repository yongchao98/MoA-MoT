import math

def solve_emperor_tomb_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total engraved characters (K).
    """
    # Define material and artifact dimensions
    material_width = 140
    material_height = 110
    square_side = 10
    circle_radius = 20
    
    # A circle with a 20cm radius requires a 40x40cm square bounding box to be cut from.
    circle_bbox_side = circle_radius * 2

    # Define character counts per artifact
    square_chars = 4
    circle_chars = 999

    # --- Step 1: Maximize the number of circles (M) ---
    # We prioritize circles due to their high character value (999 vs 4).
    # We calculate how many 40x40cm bounding boxes fit into the 140x110cm material.
    
    num_circles_along_width = material_width // circle_bbox_side
    num_circles_along_height = material_height // circle_bbox_side
    M = num_circles_along_width * num_circles_along_height

    print(f"The material is {material_width}x{material_height}cm. Each circle needs a {circle_bbox_side}x{circle_bbox_side}cm area.")
    print(f"Maximum circles that can fit: {num_circles_along_width} along the width and {num_circles_along_height} along the height.")
    print(f"Total maximum circles (M) = {num_circles_along_width} * {num_circles_along_height} = {M}")

    # --- Step 2: Calculate available squares (N) from the remaining material ---
    # The grid of circle bounding boxes creates a large used rectangle.
    used_width = num_circles_along_width * circle_bbox_side
    used_height = num_circles_along_height * circle_bbox_side

    # The remaining material forms an L-shape, which we can split into two rectangles.
    # Rectangle A: The strip remaining along the material's full width.
    rem_A_w = material_width
    rem_A_h = material_height - used_height
    squares_from_A = (rem_A_w // square_side) * (rem_A_h // square_side)

    # Rectangle B: The strip remaining next to the used area.
    rem_B_w = material_width - used_width
    rem_B_h = used_height
    squares_from_B = (rem_B_w // square_side) * (rem_B_h // square_side)

    N = squares_from_A + squares_from_B
    print(f"\nAfter cutting the circles, two rectangular areas remain:")
    print(f"1. A {rem_A_w}x{rem_A_h}cm piece, yielding {squares_from_A} squares.")
    print(f"2. A {rem_B_w}x{rem_B_h}cm piece, yielding {squares_from_B} squares.")
    print(f"Total squares from remaining area (N) = {squares_from_A} + {squares_from_B} = {N}")

    # --- Step 3: Calculate the total number of characters (K) ---
    K = (N * square_chars) + (M * circle_chars)

    # --- Step 4: Output the final result ---
    print("\n--- Final Calculation ---")
    print(f"The final equation is: K = (N * characters_per_square) + (M * characters_per_circle)")
    print(f"K = ({N} * {square_chars}) + ({M} * {circle_chars})")
    print(f"K = {N * square_chars} + {M * circle_chars}")
    print(f"K = {K}")
    
    print("\nThe optimal production ratio N:M:K is:")
    print(f"{N}:{M}:{K}")

solve_emperor_tomb_engraving()