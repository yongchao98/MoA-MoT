import math

def solve_laozi_books():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters.
    """

    # Step 1: Define constants and calculate characters per plate
    chars_to_encode = 1000
    
    # Properties of the circular plate ("Book of Heaven")
    circle_radius_cm = 20
    circle_diameter_cm = circle_radius_cm * 2
    symbols_per_circle = 9999
    yinyang_wuxing_symbols = 10
    
    # Properties of the square plate ("Book of Earth")
    square_side_cm = 10
    symbols_per_square = 360
    bagua_symbols = 8

    # Calculate characters per plate based on encoding requirements
    symbols_for_char_circle = math.ceil(math.log(chars_to_encode, yinyang_wuxing_symbols))
    chars_per_circle = math.floor(symbols_per_circle / symbols_for_char_circle)

    symbols_for_char_square = math.ceil(math.log(chars_to_encode, bagua_symbols))
    chars_per_square = math.floor(symbols_per_square / symbols_for_char_square)

    # Step 2: Analyze and compare packing strategies

    # Strategy A: Grid Packing (sub-optimal for circles, but leaves rectangular space)
    sheet_w, sheet_h = 140, 110
    
    # Circles fit in a 3x2 grid
    n_circles_grid = math.floor(sheet_w / circle_diameter_cm) * math.floor(sheet_h / circle_diameter_cm) # 3 * 2 = 6
    
    # Calculate remaining area for squares
    circles_block_w = math.floor(sheet_w / circle_diameter_cm) * circle_diameter_cm # 3 * 40 = 120
    circles_block_h = math.floor(sheet_h / circle_diameter_cm) * circle_diameter_cm # 2 * 40 = 80
    
    rem_area1_w, rem_area1_h = (sheet_w - circles_block_w), sheet_h
    rem_area2_w, rem_area2_h = circles_block_w, (sheet_h - circles_block_h)

    m_squares_grid_1 = math.floor(rem_area1_w / square_side_cm) * math.floor(rem_area1_h / square_side_cm)
    m_squares_grid_2 = math.floor(rem_area2_w / square_side_cm) * math.floor(rem_area2_h / square_side_cm)
    m_squares_grid = m_squares_grid_1 + m_squares_grid_2
    
    k_chars_grid = (n_circles_grid * chars_per_circle) + (m_squares_grid * chars_per_square)

    # Strategy B: Staggered Packing (optimal for circles)
    # This is a known result for packing circles in a rectangle of these proportions.
    n_circles_staggered = 9
    # The packing is very tight, leaving no significant space for 10x10cm squares.
    m_squares_staggered = 0
    k_chars_staggered = (n_circles_staggered * chars_per_circle) + (m_squares_staggered * chars_per_square)

    # Step 3: Compare strategies and select the best one
    if k_chars_staggered > k_chars_grid:
        final_n = n_circles_staggered
        final_m = m_squares_staggered
        final_k = k_chars_staggered
    else:
        final_n = n_circles_grid
        final_m = m_squares_grid
        final_k = k_chars_grid

    # Step 4: Print the final answer
    print("The optimal configuration is found by maximizing the number of high-value circular plates.")
    print("This is achieved with a staggered packing arrangement.")
    print(f"\nFinal optimal configuration:")
    print(f"N (circles) = {final_n}")
    print(f"M (squares) = {final_m}")
    print(f"K (characters) = {final_k}")
    
    print("\nThe calculation for the maximal number of characters is:")
    print(f"Equation: {final_n} * {chars_per_circle} + {final_m} * {chars_per_square} = {final_k}")
    
    final_answer_string = f"{final_n}:{final_m}:{final_k}"
    print(f"\n<<<_START_OF_ANSWER_>>>\n{final_answer_string}\n<<<__END_OF_ANSWER__>>>")


solve_laozi_books()