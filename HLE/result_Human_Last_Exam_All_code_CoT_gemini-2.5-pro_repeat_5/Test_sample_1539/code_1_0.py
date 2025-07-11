import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and square plates to maximize engraved characters.
    """
    # Define constants based on the problem description
    material_w = 140
    material_h = 110
    
    circle_radius = 20
    circle_diameter = 2 * circle_radius
    circle_value = 9999
    
    square_side = 10
    square_value = 360

    print("--- Plan ---")
    print("1. Prioritize cutting the high-value circular plates based on value density analysis.")
    print("2. Calculate the maximum number of circles that can fit in a grid on the 140x110cm material.")
    print("3. After allocating space for the circles, calculate the remaining area.")
    print("4. The remaining area consists of large rectangular strips and smaller spaces between the circles.")
    print("5. Calculate how many square plates can be cut from these remaining areas.")
    print("6. Sum the total number of characters from all plates to find the maximum.")
    print("\n--- Calculation ---")

    # Step 1: Maximize the number of circles in a grid layout.
    # The material is 140x110. A circle needs a 40x40cm bounding box.
    num_circles_along_w = material_w // circle_diameter
    num_circles_along_h = material_h // circle_diameter
    N = num_circles_along_w * num_circles_along_h
    
    print(f"A circle's diameter is {circle_diameter}cm. A square's side is {square_side}cm.")
    print(f"On a {material_w}x{material_h}cm sheet, we can fit {num_circles_along_w} circles along the width and {num_circles_along_h} along the height.")
    print(f"Maximum number of circular plates (N) = {num_circles_along_w} * {num_circles_along_h} = {N}")

    # Step 2: Calculate the area occupied by the circles' grid and the leftover rectangular strips.
    circles_grid_w = num_circles_along_w * circle_diameter
    circles_grid_h = num_circles_along_h * circle_diameter

    # The leftover space is an L-shape, which can be divided into two rectangles.
    # Rectangle A: The vertical part of the 'L'
    rect_A_w = material_w - circles_grid_w
    rect_A_h = material_h
    squares_from_A = (rect_A_w // square_side) * (rect_A_h // square_side)
    print(f"Remaining rectangular strip 1 has dimensions {rect_A_w}x{rect_A_h}cm, fitting {squares_from_A} squares.")

    # Rectangle B: The horizontal part of the 'L'
    rect_B_w = circles_grid_w 
    rect_B_h = material_h - circles_grid_h
    squares_from_B = (rect_B_w // square_side) * (rect_B_h // square_side)
    print(f"Remaining rectangular strip 2 has dimensions {rect_B_w}x{rect_B_h}cm, fitting {squares_from_B} squares.")
    
    # Step 3: Calculate squares from the space *between* the circles.
    # In a grid, the space between 4 circles can often fit smaller items.
    # The number of such 4-circle intersections is (num_w - 1) * (num_h - 1).
    # The space between 4 circles (diameter 40) can fit a square of side 10.
    interstitial_squares = (num_circles_along_w - 1) * (num_circles_along_h - 1)
    print(f"The {num_circles_along_w}x{num_circles_along_h} grid of circles creates {interstitial_squares} intersection(s) large enough for a 10x10cm square.")

    # Step 4: Calculate total squares (M) and total characters (K).
    M = squares_from_A + squares_from_B + interstitial_squares
    K = (N * circle_value) + (M * square_value)
    
    print("\n--- Results ---")
    print(f"Total number of circular plates (N): {N}")
    print(f"Total number of squared plates (M): {squares_from_A} + {squares_from_B} + {interstitial_squares} = {M}")
    
    print("\nFinal equation for total characters (K):")
    # The problem asks to output each number in the final equation.
    print(f"{N} * {circle_value} + {M} * {square_value} = {K}")

    # Final answer in the required format N:M:K
    final_answer = f"{N}:{M}:{K}"
    print(f"\nThe final answer is in the format N:M:K")
    print(final_answer)
    print(f"<<<{final_answer}>>>")

solve_laozi_puzzle()