import math

def solve_emperor_tomb_problem():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # Step 1: Define constants from the problem statement
    sheet_width = 140
    sheet_height = 110
    square_side = 10
    circle_radius = 20
    circle_diameter = 2 * circle_radius

    # Characters/Symbols per artifact
    chars_per_square = 4
    symbols_per_circle = 999
    unique_bio_chars = 7 * 7 * 9
    bagua_symbols_count = 8

    # Step 2: Determine the number of characters per circle
    # To uniquely represent 'unique_bio_chars' using a system of 'bagua_symbols_count' symbols,
    # we need to find n such that `8^n >= 441`. This is n >= log_8(441).
    symbols_per_char = math.ceil(math.log(unique_bio_chars, bagua_symbols_count))
    chars_per_circle = symbols_per_circle // symbols_per_char

    # Step 3: Solve the 2D packing problem
    # For cutting purposes, a circle requires its square bounding box.
    square_bounding_box = square_side
    circle_bounding_box = circle_diameter

    # Compare character density to prioritize.
    # val_density_circle = 333 / (40*40) = 0.208
    # val_density_square = 4 / (10*10) = 0.04
    # The density for circles is higher, so we prioritize maximizing M.

    # Calculate max M in a grid layout.
    m_count_in_width = sheet_width // circle_bounding_box
    m_count_in_height = sheet_height // circle_bounding_box
    M = m_count_in_width * m_count_in_height

    # This maximum number of circles is achieved with a 3x2 grid, which creates
    # a 120cm x 80cm block to be cut.
    circles_block_w = m_count_in_width * circle_bounding_box
    circles_block_h = m_count_in_height * circle_bounding_box

    # We place this block in a corner of the 140x110 sheet, leaving an L-shaped
    # area that can be divided into two rectangles.
    
    # Leftover rectangle 1
    rem_rect1_w = sheet_width - circles_block_w
    rem_rect1_h = sheet_height
    squares_in_rect1 = (rem_rect1_w // square_side) * (rem_rect1_h // square_side)

    # Leftover rectangle 2
    rem_rect2_w = circles_block_w
    rem_rect2_h = sheet_height - circles_block_h
    squares_in_rect2 = (rem_rect2_w // square_side) * (rem_rect2_h // square_side)

    N = squares_in_rect1 + squares_in_rect2

    # Step 4: Calculate the maximal number of characters K
    K = (chars_per_square * N) + (chars_per_circle * M)

    # Print the explanation and final answer
    print("To maximize the number of engraved Chinese characters, we follow this plan:")
    print("\n1. Determine Characters per Artifact:")
    print(f"- Each square for the Emperor's name holds {chars_per_square} characters.")
    print(f"- To encode {unique_bio_chars} unique bio characters using {bagua_symbols_count} symbols, we need ceil(log_8({unique_bio_chars})) = {symbols_per_char} symbols per character.")
    print(f"- Each circle with {symbols_per_circle} symbols holds {symbols_per_circle} / {symbols_per_char} = {chars_per_circle} characters.")
    
    print("\n2. Optimize Packing on the 140x110cm Sheet:")
    print("- Circles provide more characters per area, so we prioritize them.")
    print(f"- Maximum circles (M) in a grid = floor(140/{circle_bounding_box}) * floor(110/{circle_bounding_box}) = {m_count_in_width} * {m_count_in_height} = {M}.")
    print("- This uses a 120x80cm area. The remaining L-shaped space is used for squares.")
    print("- The leftover space yields two rectangular areas for squares:")
    print(f"  - Area 1 (20x110cm): floor(20/{square_side}) * floor(110/{square_side}) = {squares_in_rect1} squares.")
    print(f"  - Area 2 (120x30cm): floor(120/{square_side}) * floor(30/{square_side}) = {squares_in_rect2} squares.")
    print(f"- Total squares (N) = {squares_in_rect1} + {squares_in_rect2} = {N}.")

    print("\n3. Calculate Final Result:")
    print("The maximum total characters (K) is calculated as:")
    print(f"K = (N × Characters per Square) + (M × Characters per Circle)")
    print(f"K = ({N} × {chars_per_square}) + ({M} × {chars_per_circle}) = {K}")
    
    print(f"\nFinal Answer (N:M:K) = {N}:{M}:{K}")

solve_emperor_tomb_problem()
<<<58:6:2230>>>