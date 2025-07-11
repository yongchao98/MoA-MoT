import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters from a rectangular meteorite sheet.
    """
    # --- Define Constants ---
    rect_width = 140  # cm
    rect_height = 110 # cm

    circle_radius = 20 # cm
    circle_diameter = circle_radius * 2
    circle_chars = 9999

    square_side = 10 # cm
    square_chars = 360

    # --- Strategy: Maximize the number of high-value circular plates first ---
    # We place the circles in a grid. We check both orientations of the rectangle
    # to see which one yields a better result, though in this case they are equivalent.

    # Orientation 1: 140x110
    circles_along_width = rect_width // circle_diameter
    circles_along_height = rect_height // circle_diameter
    num_circles_1 = circles_along_width * circles_along_height

    # Calculate remaining space for Orientation 1
    used_width_1 = circles_along_width * circle_diameter
    used_height_1 = circles_along_height * circle_diameter

    # The remaining area consists of two rectangular strips
    rem_strip_A_w = rect_width
    rem_strip_A_h = rect_height - used_height_1
    rem_strip_B_w = rect_width - used_width_1
    rem_strip_B_h = used_height_1

    squares_in_A = (rem_strip_A_w // square_side) * (rem_strip_A_h // square_side)
    squares_in_B = (rem_strip_B_w // square_side) * (rem_strip_B_h // square_side)
    num_squares_1 = squares_in_A + squares_in_B

    total_chars_1 = (num_circles_1 * circle_chars) + (num_squares_1 * square_chars)

    # For this problem, rotating the sheet (110x140) yields the same number of circles
    # and the same total area, which also results in the same number of squares.
    # So we can proceed with the calculated values.
    N = num_circles_1
    M = num_squares_1
    K = total_chars_1

    # --- Output the results ---
    print("To maximize the number of characters, we should prioritize the circular plates.")
    print(f"The 140x110cm material can fit a grid of {circles_along_width}x{circles_along_height} circular plates.")
    print(f"Number of circular plates (N): {N}")
    print(f"This leaves space for {M} square plates.")
    print(f"Number of square plates (M): {M}")
    print("\nThe final calculation for the maximum number of characters (K) is:")
    print(f"{N} * {circle_chars} + {M} * {square_chars} = {K}")
    print("\nThe final answer in the format N:M:K is:")
    print(f"{N}:{M}:{K}")

solve_laozi_engraving()
<<<6:58:80874>>>