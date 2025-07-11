import math

def solve_engraving_problem():
    """
    Calculates the maximum number of Chinese characters that can be engraved
    by optimizing the cutting of squares and circles from a rectangular sheet.
    """
    # Material and item dimensions
    sheet_width = 140
    sheet_height = 110
    square_side = 10
    circle_diameter = 40  # 20cm radius

    # Characters per item
    chars_per_square = 4
    chars_per_circle = 999

    # --- Step 1: Maximize the number of circles (M) ---
    # The contribution of circles (999 chars) is far greater than squares (4 chars),
    # so we prioritize fitting as many circles as possible.

    # Calculate how many circles fit along each dimension of the sheet
    circles_along_width = sheet_width // circle_diameter
    circles_along_height = sheet_height // circle_diameter
    
    # The maximum number of circles is the product of these counts
    M = circles_along_width * circles_along_height

    # --- Step 2: Maximize squares (N) in the remaining area ---
    # The circles occupy a rectangular block on the sheet
    circles_occupied_width = circles_along_width * circle_diameter
    circles_occupied_height = circles_along_height * circle_diameter

    # The remaining space can be divided into two rectangular strips
    # Strip 1:
    rem_strip1_width = sheet_width - circles_occupied_width
    rem_strip1_height = sheet_height 
    squares_in_strip1 = (rem_strip1_width // square_side) * (rem_strip1_height // square_side)
    
    # Strip 2:
    rem_strip2_width = circles_occupied_width
    rem_strip2_height = sheet_height - circles_occupied_height
    squares_in_strip2 = (rem_strip2_width // square_side) * (rem_strip2_height // square_side)

    # Total number of squares
    N = squares_in_strip1 + squares_in_strip2

    # --- Step 3: Calculate the maximum total characters (K) ---
    K = (chars_per_square * N) + (chars_per_circle * M)

    # --- Step 4: Display the results ---
    print("Optimization Plan:")
    print(f"- First, maximize the number of circles (M) on the {sheet_width}x{sheet_height}cm sheet.")
    print(f"- Each circle needs a {circle_diameter}x{circle_diameter}cm space.")
    print(f"- We can fit {circles_along_width} along the width and {circles_along_height} along the height, for a total of M = {M} circles.")
    print("- Then, we fit squares in the remaining area.")
    print(f"- After placing the circles, the remaining area yields N = {N} squares.")
    print("\nFinal Calculation:")
    print(f"The equation for the total characters is K = (characters per square) * N + (characters per circle) * M")
    print(f"K = {chars_per_square} * {N} + {chars_per_circle} * {M} = {K}")
    
    print("\nResult:")
    print(f"The optimal arrangement is N:M:K = {N}:{M}:{K}")

solve_engraving_problem()