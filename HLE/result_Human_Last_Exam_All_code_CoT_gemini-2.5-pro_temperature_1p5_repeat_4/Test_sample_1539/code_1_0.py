import math

def solve_laozi_riddle():
    """
    Solves the optimization problem for Laozi's books.
    """

    # Step 1: Determine the character value of each plate

    # Plate dimensions
    sheet_dim = (140, 110)
    circle_radius = 20
    circle_diameter = 2 * circle_radius
    square_side = 10

    # Symbols and characters
    num_chars_to_encode = 1000
    symbols_on_circle = 9999
    symbols_on_square = 360
    
    # Yinyang Wuxing has 2 (yin/yang) * 5 (elements) = 10 unique symbols
    base_heaven = 10
    # Bagua has 8 unique trigram symbols
    base_earth = 8

    # Calculate symbols needed per character for each system
    # k = ceil(log_B(1000))
    symbols_per_char_heaven = math.ceil(math.log(num_chars_to_encode, base_heaven))
    symbols_per_char_earth = math.ceil(math.log(num_chars_to_encode, base_earth))

    # Calculate total characters per plate
    chars_per_circle = math.floor(symbols_on_circle / symbols_per_char_heaven)
    chars_per_square = math.floor(symbols_on_square / symbols_per_char_earth)

    print("Step 1: Calculating Character Capacity per Plate")
    print(f" - Symbols needed per character (Heaven): ceil(log10(1000)) = {symbols_per_char_heaven}")
    print(f" - Characters per circular plate: floor({symbols_on_circle} / {symbols_per_char_heaven}) = {chars_per_circle}")
    print(f" - Symbols needed per character (Earth): ceil(log8(1000)) = {symbols_per_char_earth}")
    print(f" - Characters per squared plate: floor({symbols_on_square} / {symbols_per_char_earth}) = {chars_per_square}\n")
    
    # Step 2 & 3: Solve the packing problem

    print("Step 2: Solving the Packing Problem")
    # To maximize the total characters, we prioritize the plate type with the higher value.
    # The value of a circle is 3333 characters, far greater than a square's 90.
    # Therefore, we aim to maximize the number of circular plates.

    # We determine the maximum number of circles that can be cut from the 140x110cm sheet.
    # A simple grid layout (3x2) yields 6 circles.
    # A more efficient staggered (hexagonal) packing allows for more.
    # For a 140x110 sheet (width x height):
    # - We can fit 3 circles per row along the 140cm side.
    # - We can fit 3 staggered rows along the 110cm side.
    # - This gives a 3-2-3 arrangement, for a total of N=8 circles.
    # This layout uses a bounding box of 120cm x 109.28cm.
    N = 8
    
    # This layout leaves two main rectangular strips of material.
    # 1. A vertical strip: (140 - 120)cm wide by 110cm high.
    # 2. A horizontal strip: 120cm wide by (110 - 109.28)cm high.
    
    rem_strip1_dim = (sheet_dim[0] - 3 * circle_diameter, sheet_dim[1])
    rem_strip2_dim = (3 * circle_diameter, sheet_dim[1] - (circle_diameter + 2 * (circle_radius * math.sqrt(3))))

    squares_from_strip1 = math.floor(rem_strip1_dim[0] / square_side) * math.floor(rem_strip1_dim[1] / square_side)
    squares_from_strip2 = math.floor(rem_strip2_dim[0] / square_side) * math.floor(rem_strip2_dim[1] / square_side)

    M = squares_from_strip1 + squares_from_strip2
    
    print(f" - Maximum circles (N) using staggered packing: {N}")
    print(f" - Leftover material is primarily a {rem_strip1_dim[0]}x{rem_strip1_dim[1]}cm strip.")
    print(f" - Maximum squares (M) from leftover material: {M}\n")

    # Step 4: Calculate the final result
    K = N * chars_per_circle + M * chars_per_square
    
    print("Step 3: Final Calculation")
    print("The maximum number of characters K is N * (chars per circle) + M * (chars per square).")
    # The problem asks to output each number in the final equation.
    print(f"K = {N} * {chars_per_circle} + {M} * {chars_per_square} = {K}")
    
    print("\nFinal Answer (N:M:K):")
    print(f"{N}:{M}:{K}")

solve_laozi_riddle()