# Define constants for the problem
MATERIAL_WIDTH = 140
MATERIAL_HEIGHT = 110
SQUARE_SIDE = 10
CIRCLE_DIAMETER = 40  # A circle of 20cm radius fits in a 40x40cm square
CHARS_PER_SQUARE = 4
CHARS_PER_CIRCLE = 999

def solve_engraving_problem():
    """
    Calculates the optimal number of squares and circles to maximize engraved characters.
    """
    print("Step 1: Determine the optimal strategy.")
    # Calculate characters per 40x40cm area for both options
    squares_in_circle_area = (CIRCLE_DIAMETER // SQUARE_SIDE) ** 2
    chars_from_squares_in_area = squares_in_circle_area * CHARS_PER_SQUARE
    print(f"A 40x40 cm area can fit 1 circle for {CHARS_PER_CIRCLE} characters.")
    print(f"The same area can fit {squares_in_circle_area} squares for {chars_from_squares_in_area} characters.")
    print("Strategy: Prioritize circles due to much higher character count.\n")

    print("Step 2: Calculate the maximum number of circles (M).")
    # We are packing 40x40cm squares into a 140x110cm rectangle.
    num_circles_along_width = MATERIAL_WIDTH // CIRCLE_DIAMETER
    num_circles_along_height = MATERIAL_HEIGHT // CIRCLE_DIAMETER
    M = num_circles_along_width * num_circles_along_height
    print(f"Circles along 140cm side: floor(140 / 40) = {num_circles_along_width}")
    print(f"Circles along 110cm side: floor(110 / 40) = {num_circles_along_height}")
    print(f"Maximum number of circles (M) = {num_circles_along_width} * {num_circles_along_height} = {M}\n")

    print("Step 3: Calculate the number of squares (N) in the remaining area.")
    # The 6 circles are arranged in a 3x2 grid, occupying a 120x80cm area.
    used_width = num_circles_along_width * CIRCLE_DIAMETER
    used_height = num_circles_along_height * CIRCLE_DIAMETER
    print(f"Area used by circles: {used_width}cm x {used_height}cm.")

    # Calculate the remaining area. This forms an L-shape.
    # We can split this L-shape into two rectangles to calculate the number of squares.
    # Rectangle 1: (MATERIAL_WIDTH) x (MATERIAL_HEIGHT - used_height)
    # Rectangle 2: (used_height) x (MATERIAL_WIDTH - used_width) - typo, it's (used_width) x ...
    # Corrected: It's the remaining L-shape.
    # Rect 1: The strip along the full width: 140 x (110-80) = 140 x 30
    # Rect 2: The strip along the side of the used area: (140-120) x 80 = 20 x 80
    
    # Let's use area calculation, which is simpler as the squares fit perfectly.
    total_area = MATERIAL_WIDTH * MATERIAL_HEIGHT
    used_area = used_width * used_height
    remaining_area = total_area - used_area
    square_area = SQUARE_SIDE ** 2
    N = remaining_area // square_area
    
    print(f"Total material area: {total_area} cm^2.")
    print(f"Area used by circles: {used_area} cm^2.")
    print(f"Remaining area for squares: {remaining_area} cm^2.")
    print(f"Number of 10x10 squares (N) = {remaining_area} / {square_area} = {N}\n")

    print("Step 4: Calculate the total number of characters (K).")
    K = (CHARS_PER_SQUARE * N) + (CHARS_PER_CIRCLE * M)
    print(f"The formula is K = ({CHARS_PER_SQUARE} * N) + ({CHARS_PER_CIRCLE} * M)")
    print(f"Calculation: K = ({CHARS_PER_SQUARE} * {N}) + ({CHARS_PER_CIRCLE} * {M}) = {CHARS_PER_SQUARE * N} + {CHARS_PER_CIRCLE * M} = {K}\n")
    
    print("--- Final Result ---")
    print(f"The optimal combination is N={N} squares and M={M} circles.")
    print(f"This yields a maximum of K={K} characters.")
    print(f"Answer in N:M:K format: {N}:{M}:{K}")

solve_engraving_problem()