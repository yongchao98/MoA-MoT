import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters.
    """
    # --- Problem Parameters ---
    sheet_width = 140
    sheet_height = 110
    
    # Circle properties
    circle_radius = 20
    circle_chars = 9999
    
    # Square properties
    square_side = 10
    square_chars = 360

    # --- Step 1: Justify the Greedy Strategy ---
    # The problem is to maximize K = N * 9999 + M * 360.
    # A circle requires a 40x40cm bounding box to be cut.
    # A square requires a 10x10cm box.
    # One 40x40 box could instead be used to make (40*40)/(10*10) = 16 squares.
    # Value of 1 circle = 9999 characters.
    # Value of 16 squares = 16 * 360 = 5760 characters.
    # Since 9999 > 5760, it's always better to make a circle.
    # Therefore, we will maximize the number of circles first.

    # --- Step 2: Calculate the maximum number of circles (N) ---
    circle_diameter = 2 * circle_radius
    
    # Calculate how many circles fit along the width and height
    num_circles_along_width = sheet_width // circle_diameter
    num_circles_along_height = sheet_height // circle_diameter
    
    N = num_circles_along_width * num_circles_along_height

    # --- Step 3: Calculate the number of squares (M) from the remaining area ---
    # Calculate the dimensions of the rectangular area used by the circles
    width_used_by_circles = num_circles_along_width * circle_diameter
    height_used_by_circles = num_circles_along_height * circle_diameter
    
    # Calculate the area of the entire sheet and the area used by circles
    total_sheet_area = sheet_width * sheet_height
    area_used_by_circles = width_used_by_circles * height_used_by_circles
    
    # The remaining area can be perfectly tiled by squares
    remaining_area = total_sheet_area - area_used_by_circles
    square_area = square_side * square_side
    M = remaining_area // square_area

    # --- Step 4: Calculate the total number of characters (K) ---
    K = (N * circle_chars) + (M * square_chars)

    # --- Step 5: Output the results ---
    print(f"Based on the 140x110cm material:")
    print(f"The optimal number of circular plates (N) is: {N}")
    print(f"The optimal number of squared plates (M) is: {M}")
    print(f"The maximal number of engraved characters (K) is: {K}")
    print("\nThis is calculated as follows:")
    print(f"({N} circular plates * {circle_chars} chars/plate) + ({M} squared plates * {square_chars} chars/plate) = {K} characters")
    
    print("\nFinal Equation:")
    print(f"{N} * {circle_chars} + {M} * {square_chars} = {K}")
    
    print("\nFinal Answer (N:M:K):")
    print(f"{N}:{M}:{K}")

solve_laozi_puzzle()
<<<6:58:80874>>>