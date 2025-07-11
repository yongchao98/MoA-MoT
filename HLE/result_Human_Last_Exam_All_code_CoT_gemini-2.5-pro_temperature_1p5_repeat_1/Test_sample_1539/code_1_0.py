import math

def solve_engraving_maximization():
    """
    Calculates the optimal number of circular and square plates to maximize
    the total number of engraved characters from a rectangular sheet of material.
    """
    # Define the parameters from the problem description
    sheet_width = 140
    sheet_height = 110

    # Parameters for the "Book of Heaven" circular plates
    circle_radius = 20
    # The space required for a circle is its bounding box (a square with side = diameter)
    circle_bounding_box_side = circle_radius * 2
    circle_chars = 9999

    # Parameters for the "Book of Earth" square plates
    square_side = 10
    square_chars = 360

    # Step 1: Prioritize packing the most "valuable" plates.
    # A circular plate uses a 40x40=1600 cm^2 area for 9999 chars. (Value density: ~6.25 chars/cm^2)
    # 16 square plates would use the same 1600 cm^2 area for 16*360=5760 chars.
    # Therefore, circular plates are more valuable and should be maximized first.

    # Step 2: Calculate the maximum number of circular plates (N).
    # We fit the 40x40 cm bounding boxes into the 140x110 cm sheet.
    circles_along_width = sheet_width // circle_bounding_box_side
    circles_along_height = sheet_height // circle_bounding_box_side
    
    # Total number of circular plates
    N = circles_along_width * circles_along_height

    # Step 3: Calculate remaining area and pack square plates (M).
    # The packed circular plates form a rectangular block.
    circles_block_width = circles_along_width * circle_bounding_box_side
    circles_block_height = circles_along_height * circle_bounding_box_side

    # This leaves an L-shaped area, which can be split into two rectangles.
    # Rectangle 1
    rem_rect1_width = sheet_width - circles_block_width
    rem_rect1_height = sheet_height
    squares_in_rect1 = (rem_rect1_width // square_side) * (rem_rect1_height // square_side)

    # Rectangle 2
    rem_rect2_width = circles_block_width
    rem_rect2_height = sheet_height - circles_block_height
    squares_in_rect2 = (rem_rect2_width // square_side) * (rem_rect2_height // square_side)

    # Total number of square plates
    M = squares_in_rect1 + squares_in_rect2

    # Step 4: Calculate the maximal number of characters (K).
    total_chars_from_circles = N * circle_chars
    total_chars_from_squares = M * square_chars
    K = total_chars_from_circles + total_chars_from_squares

    # Print the explanation and the final result
    print("To maximize the engraved characters, the optimal combination is:")
    print(f"N = {N} circular plates")
    print(f"M = {M} square plates")
    print("\nThe calculation for the total characters (K) is as follows:")
    print(f"K = (N * chars_per_circle) + (M * chars_per_square)")
    print(f"K = ({N} * {circle_chars}) + ({M} * {square_chars})")
    print(f"K = {total_chars_from_circles} + {total_chars_from_squares}")
    print(f"K = {K}")
    print("\nSo, the final answer in the N:M:K format is:")
    print(f"<<<{N}:{M}:{K}>>>")

# Run the solver
solve_engraving_maximization()