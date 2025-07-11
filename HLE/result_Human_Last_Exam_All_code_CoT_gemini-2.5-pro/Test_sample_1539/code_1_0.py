import math

def solve_laozi_puzzle():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved Chinese characters.
    """
    # Step 1 & 2: Define constants and calculate characters per plate
    
    # Material and plate dimensions
    MATERIAL_W = 140  # cm
    MATERIAL_H = 110  # cm
    CIRCLE_RADIUS = 20  # cm
    SQUARE_SIDE = 10   # cm

    # A circle with radius 20cm needs a 40x40cm square area to be cut out
    CIRCLE_CUT_SIDE = CIRCLE_RADIUS * 2

    # Symbol and character information
    SYMBOLS_HEAVEN = 9999 # yinyang wuxing symbols per circular plate
    SYMBOLS_EARTH = 360   # bagua symbols per square plate
    NUM_CHINESE_CHARS = 1000

    # There are 10 unique yinyang wuxing symbols (2 yin/yang * 5 wuxing)
    SYMBOLS_PER_CHAR_HEAVEN = math.ceil(math.log(NUM_CHINESE_CHARS, 10))
    # There are 8 unique bagua symbols
    SYMBOLS_PER_CHAR_EARTH = math.ceil(math.log(NUM_CHINESE_CHARS, 8))

    # Calculate how many characters fit on one plate
    CHARS_PER_HEAVEN_PLATE = math.floor(SYMBOLS_HEAVEN / SYMBOLS_PER_CHAR_HEAVEN)
    CHARS_PER_EARTH_PLATE = math.floor(SYMBOLS_EARTH / SYMBOLS_PER_CHAR_EARTH)

    # Step 3: Solve the packing problem
    # We prioritize fitting the more valuable circular plates first.
    # We will use a simple grid packing strategy.
    
    # Check orientation 1: 140x110
    n_circles_w1 = MATERIAL_W // CIRCLE_CUT_SIDE
    n_circles_h1 = MATERIAL_H // CIRCLE_CUT_SIDE
    n_circles1 = n_circles_w1 * n_circles_h1
    
    # Calculate remaining area for squares
    used_w1 = n_circles_w1 * CIRCLE_CUT_SIDE
    used_h1 = n_circles_h1 * CIRCLE_CUT_SIDE
    
    # Area 1: The strip left on the side
    rem_w1 = MATERIAL_W - used_w1
    squares_in_rem1 = (rem_w1 // SQUARE_SIDE) * (MATERIAL_H // SQUARE_SIDE)
    
    # Area 2: The strip left on the bottom (on the area already used in width)
    rem_h2 = MATERIAL_H - used_h1
    squares_in_rem2 = (used_w1 // SQUARE_SIDE) * (rem_h2 // SQUARE_SIDE)
    
    n_squares1 = squares_in_rem1 + squares_in_rem2
    
    total_chars1 = n_circles1 * CHARS_PER_HEAVEN_PLATE + n_squares1 * CHARS_PER_EARTH_PLATE

    # Check orientation 2: 110x140 (swap width and height)
    n_circles_w2 = MATERIAL_H // CIRCLE_CUT_SIDE
    n_circles_h2 = MATERIAL_W // CIRCLE_CUT_SIDE
    n_circles2 = n_circles_w2 * n_circles_h2

    used_w2 = n_circles_w2 * CIRCLE_CUT_SIDE
    used_h2 = n_circles_h2 * CIRCLE_CUT_SIDE
    
    rem_w1_orient2 = MATERIAL_H - used_w2
    squares_in_rem1_orient2 = (rem_w1_orient2 // SQUARE_SIDE) * (MATERIAL_W // SQUARE_SIDE)
    
    rem_h2_orient2 = MATERIAL_W - used_h2
    squares_in_rem2_orient2 = (used_w2 // SQUARE_SIDE) * (rem_h2_orient2 // SQUARE_SIDE)
    
    n_squares2 = squares_in_rem1_orient2 + squares_in_rem2_orient2
    
    total_chars2 = n_circles2 * CHARS_PER_HEAVEN_PLATE + n_squares2 * CHARS_PER_EARTH_PLATE

    # Find the best result
    if total_chars1 > total_chars2:
        N = n_circles1
        M = n_squares1
        K = total_chars1
    else:
        # Note: In this specific case, both orientations yield the same result.
        N = n_circles2
        M = n_squares2
        K = total_chars2

    # Step 4: Print the final result
    print(f"Optimal number of circular plates (N): {N}")
    print(f"Optimal number of squared plates (M): {M}")
    print(f"Maximum number of characters (K): {K}")
    print("\nCalculation:")
    print(f"{N} * {CHARS_PER_HEAVEN_PLATE} + {M} * {CHARS_PER_EARTH_PLATE} = {K}")
    print("\nFinal Answer in N:M:K format:")
    print(f"{N}:{M}:{K}")

solve_laozi_puzzle()
<<<6:58:25218>>>