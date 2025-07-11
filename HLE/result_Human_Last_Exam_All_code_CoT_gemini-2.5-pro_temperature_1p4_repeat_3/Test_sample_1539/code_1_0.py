import math

def solve_engraving_problem():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters from a given rectangular material.
    """
    # Define material and plate properties
    material_width = 140
    material_height = 110
    
    circle_radius = 20
    circle_diameter = circle_radius * 2
    chars_per_circle = 9999
    
    square_side = 10
    chars_per_square = 360

    print("Step 1: Determine the optimal cutting strategy.")
    # A 40x40 cm space can fit 1 circle or 16 squares.
    chars_from_circle_area = chars_per_circle
    squares_in_circle_area = (circle_diameter // square_side) * (circle_diameter // square_side)
    chars_from_square_area = squares_in_circle_area * chars_per_square
    print(f"A {circle_diameter}x{circle_diameter} cm area yields {chars_from_circle_area} characters if used for a circle.")
    print(f"The same area yields {squares_in_circle_area} * {chars_per_square} = {chars_from_square_area} characters if used for squares.")
    print("Since {} > {}, it is always better to cut a circle. Our strategy is to maximize circles first.\n".format(chars_from_circle_area, chars_from_square_area))

    print("Step 2: Calculate the maximum number of circular plates (N).")
    # Pack circles in a grid layout
    circles_along_width = material_width // circle_diameter
    circles_along_height = material_height // circle_diameter
    N = circles_along_width * circles_along_height
    print(f"Plates fit along width: {material_width} // {circle_diameter} = {circles_along_width}")
    print(f"Plates fit along height: {material_height} // {circle_diameter} = {circles_along_height}")
    print(f"Total circular plates (N) = {circles_along_width} * {circles_along_height} = {N}\n")

    # Calculate the area used by the circles' bounding boxes
    used_width = circles_along_width * circle_diameter
    used_height = circles_along_height * circle_diameter
    
    print("Step 3: Calculate the maximum number of square plates (M) from leftovers.")
    # The leftover area can be divided into two rectangles.
    # Method 1: (140-120)x110 and 120x(110-80)
    rem_rect1_w = material_width - used_width
    rem_rect1_h = material_height
    
    rem_rect2_w = used_width
    rem_rect2_h = material_height - used_height

    squares_from_rect1 = (rem_rect1_w // square_side) * (rem_rect1_h // square_side)
    squares_from_rect2 = (rem_rect2_w // square_side) * (rem_rect2_h // square_side)
    M = squares_from_rect1 + squares_from_rect2
    
    print(f"The space for circles is {used_width}x{used_height} cm.")
    print(f"This leaves two rectangular areas: {rem_rect1_w}x{rem_rect1_h} cm and {rem_rect2_w}x{rem_rect2_h} cm.")
    print(f"Squares from first rectangle: ({rem_rect1_w}//{square_side}) * ({rem_rect1_h}//{square_side}) = {squares_from_rect1}")
    print(f"Squares from second rectangle: ({rem_rect2_w}//{square_side}) * ({rem_rect2_h}//{square_side}) = {squares_from_rect2}")
    print(f"Total square plates (M) = {squares_from_rect1} + {squares_from_rect2} = {M}\n")

    print("Step 4: Calculate the total number of characters (K).")
    # Calculate total characters K
    total_chars_from_circles = N * chars_per_circle
    total_chars_from_squares = M * chars_per_square
    K = total_chars_from_circles + total_chars_from_squares
    
    print(f"The final equation is: (N * Chars/Circle) + (M * Chars/Square) = K")
    print(f"({N} * {chars_per_circle}) + ({M} * {chars_per_square}) = {total_chars_from_circles} + {total_chars_from_squares} = {K}\n")

    print("Step 5: Final Answer.")
    print(f"The optimal combination is N:M:K")
    print(f"{N}:{M}:{K}")

solve_engraving_problem()