import math

def solve_engraving_problem():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters from a rectangular sheet of material.
    """
    # Define constants based on the problem description
    sheet_length = 140
    sheet_width = 110
    
    circle_radius = 20
    circle_diameter = 2 * circle_radius
    circle_chars = 9999
    
    square_side = 10
    square_chars = 360

    # Step 1 & 2: Strategy is to maximize circles. We calculate the max number
    # of circles that fit in a simple grid layout.
    # The orientation 140x110 allows floor(140/40) x floor(110/40) = 3x2 circles.
    # The orientation 110x140 allows floor(110/40) x floor(140/40) = 2x3 circles.
    # Both result in 6 circles, so we proceed with one orientation.
    
    num_circles_along_length = sheet_length // circle_diameter
    num_circles_along_width = sheet_width // circle_diameter
    
    N = num_circles_along_length * num_circles_along_width

    # Step 3: Calculate the area occupied by the circles' bounding boxes
    circles_area_length = num_circles_along_length * circle_diameter
    circles_area_width = num_circles_along_width * circle_diameter
    
    # Step 4: Calculate how many squares can fit in the remaining L-shaped area.
    # The L-shape is split into two non-overlapping rectangles.
    
    # Rectangle 1: (full sheet length) x (remaining width)
    rem_rect1_l = sheet_length
    rem_rect1_w = sheet_width - circles_area_width
    squares_in_rem1 = (rem_rect1_l // square_side) * (rem_rect1_w // square_side)
    
    # Rectangle 2: (remaining length) x (width used by circles)
    rem_rect2_l = sheet_length - circles_area_length
    rem_rect2_w = circles_area_width
    squares_in_rem2 = (rem_rect2_l // square_side) * (rem_rect2_w // square_side)
    
    M = squares_in_rem1 + squares_in_rem2

    # Step 5: Calculate the maximum number of characters
    K = (N * circle_chars) + (M * square_chars)

    # Print the final results as requested
    print("To maximize the number of engraved characters, Laozi should produce:")
    print(f"- {N} circular plates for the 'Book of Heaven'")
    print(f"- {M} square plates for the 'Book of Earth'")
    print("\nThis combination yields the maximal number of characters.")
    print("\nThe final equation is:")
    print(f"{N} * {circle_chars} + {M} * {square_chars} = {K}")
    
    print(f"\nFinal Answer (N:M:K):")
    # Do not remove the next line, it is for the final answer.
    print(f"<<<{N}:{M}:{K}>>>")

solve_engraving_problem()