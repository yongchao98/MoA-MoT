def solve_laozi_engraving():
    """
    Calculates the maximum number of characters that can be engraved by optimizing the cutting of circular and square plates from a rectangular material.
    """

    # --- Problem Parameters ---
    # Material dimensions in cm
    rect_length = 140
    rect_width = 110

    # Plate dimensions in cm
    circle_radius = 20
    circle_diameter = circle_radius * 2
    square_side = 10

    # Character capacity per plate
    chars_per_circle = 9999
    chars_per_square = 360

    # --- Strategy: Maximize circles first due to higher character value per area ---

    # Step 1: Calculate the maximum number of circular plates (N)
    # The number of circles that fit is determined by how many of their 40x40cm bounding boxes fit.
    num_circles_along_length = rect_length // circle_diameter
    num_circles_along_width = rect_width // circle_diameter
    num_circles = num_circles_along_length * num_circles_along_width
    N = num_circles

    # Step 2: Calculate the remaining area and the maximum number of square plates (M)
    # Area occupied by the bounding boxes of the circles
    circles_occupied_length = num_circles_along_length * circle_diameter
    circles_occupied_width = num_circles_along_width * circle_diameter

    # The remaining material is calculated as two non-overlapping rectangular strips.
    # Strip 1 (the vertical remainder)
    strip1_length = rect_length - circles_occupied_length
    strip1_width = rect_width
    squares_in_strip1 = (strip1_length // square_side) * (strip1_width // square_side)

    # Strip 2 (the horizontal remainder)
    strip2_length = circles_occupied_length
    strip2_width = rect_width - circles_occupied_width
    squares_in_strip2 = (strip2_length // square_side) * (strip2_width // square_side)
    
    # Total number of squares
    num_squares = squares_in_strip1 + squares_in_strip2
    M = num_squares

    # Step 3: Calculate the total maximum number of characters (K)
    total_chars_from_circles = N * chars_per_circle
    total_chars_from_squares = M * chars_per_square
    K = total_chars_from_circles + total_chars_from_squares

    # --- Output the results ---
    print(f"The optimal number of circular plates (N) is: {N}")
    print(f"The optimal number of square plates (M) is: {M}")
    print("\nThe calculation for the maximal number of characters (K) is:")
    print(f"({N} circular plates * {chars_per_circle} chars/plate) + ({M} square plates * {chars_per_square} chars/plate)")
    print(f"= {total_chars_from_circles} + {total_chars_from_squares}")
    print(f"= {K}")
    print("\nThe final answer in the format N:M:K is:")
    print(f"<<<{N}:{M}:{K}>>>")

solve_laozi_engraving()