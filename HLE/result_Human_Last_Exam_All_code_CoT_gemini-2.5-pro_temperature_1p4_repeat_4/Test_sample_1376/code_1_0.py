import math

def solve_emperor_tomb_engraving():
    """
    This script calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K) for Emperor Qin's tomb.
    """
    # --- Problem Constants ---
    CONTAINER_W = 140
    CONTAINER_H = 110
    SQUARE_SIDE = 10
    CIRCLE_RADIUS = 20
    CIRCLE_DIAMETER = CIRCLE_RADIUS * 2
    CHARS_PER_SQUARE = 4
    SYMBOLS_PER_CIRCLE = 999
    
    # --- Step 1: Analyze Character Encoding ---
    # The bio has a set of 7x7x9 unique characters.
    unique_bio_chars = 7 * 7 * 9
    # These are encoded using Bagua symbols (base 8). We find how many symbols
    # are needed to represent each unique character.
    # We need to find the smallest integer 's' such that 8^s >= unique_bio_chars.
    symbols_per_char = math.ceil(math.log(unique_bio_chars) / math.log(8))
    # Calculate how many full characters fit on one circle.
    chars_per_circle = math.floor(SYMBOLS_PER_CIRCLE / symbols_per_char)

    print("Step 1: Analyzing the character encoding for the bio.")
    print(f"The bio has {unique_bio_chars} unique characters to encode.")
    print(f"Using 8 Bagua symbols, we need {symbols_per_char} symbols per character (since 8^{symbols_per_char} >= {unique_bio_chars}).")
    print(f"One circle with {SYMBOLS_PER_CIRCLE} symbols can hold floor({SYMBOLS_PER_CIRCLE} / {symbols_per_char}) = {chars_per_circle} bio characters.")
    print("-" * 50)

    # --- Step 2: Maximize the Number of Circles (M) ---
    # The objective is to maximize K = (4 * N) + (333 * M).
    # Since the coefficient for M is much larger, we prioritize maximizing M.
    # We pack the bounding box of the circle (40x40cm).
    circles_w = math.floor(CONTAINER_W / CIRCLE_DIAMETER)
    circles_h = math.floor(CONTAINER_H / CIRCLE_DIAMETER)
    M = circles_w * circles_h
    
    print("Step 2: Maximizing the number of circles (M).")
    print("To maximize total characters, we must first cut the maximum number of circles.")
    print(f"A 140x110cm area can fit a grid of {circles_w}x{circles_h} circles (40cm diameter).")
    print(f"Maximum number of circles (M) = {M}")
    print("-" * 50)
    
    # --- Step 3: Calculate Squares (N) in Remaining Area ---
    # The grid of circles occupies a (3*40) x (2*40) = 120x80cm area.
    packed_w = circles_w * CIRCLE_DIAMETER
    packed_h = circles_h * CIRCLE_DIAMETER
    
    # The remaining L-shaped area is split into two rectangles.
    # Rectangle 1
    rem_rect1_w = CONTAINER_W - packed_w
    rem_rect1_h = CONTAINER_H
    squares_in_rect1 = math.floor(rem_rect1_w / SQUARE_SIDE) * math.floor(rem_rect1_h / SQUARE_SIDE)
    
    # Rectangle 2
    rem_rect2_w = packed_w
    rem_rect2_h = CONTAINER_H - packed_h
    squares_in_rect2 = math.floor(rem_rect2_w / SQUARE_SIDE) * math.floor(rem_rect2_h / SQUARE_SIDE)

    N = squares_in_rect1 + squares_in_rect2

    print("Step 3: Calculating the number of squares (N) in the remaining area.")
    print(f"After cutting {M} circles in a {packed_w}x{packed_h}cm block, we fit squares in the rest.")
    print(f"The remaining space fits a maximum of {N} squares (10x10cm).")
    print("-" * 50)
    
    # --- Step 4: Final Calculation ---
    k_from_squares = CHARS_PER_SQUARE * N
    k_from_circles = chars_per_circle * M
    K = k_from_squares + k_from_circles
    
    print("Step 4: Calculating the maximum total characters (K).")
    print("The optimal production is:")
    print(f"Number of squares (N) = {N}")
    print(f"Number of circles (M) = {M}")
    print("\nFinal Result:")
    print("The maximal number of characters (K) is K = (chars from squares) + (chars from circles)")
    print(f"K = ({CHARS_PER_SQUARE} * {N}) + ({chars_per_circle} * {M})")
    print(f"K = {k_from_squares} + {k_from_circles}")
    print(f"K = {K}")
    print("\nAnswer in the format N:M:K")
    print(f"{N}:{M}:{K}")

solve_emperor_tomb_engraving()