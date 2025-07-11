import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the number of engraved Chinese characters from a given meteorite sheet.
    """
    # --- Step 1: Define constants and calculate character value per plate ---
    RECT_W, RECT_H = 140, 110
    CIRCLE_R, SQUARE_S = 20, 10
    SYMBOLS_PER_CIRCLE, SYMBOLS_PER_SQUARE = 9999, 360
    TOTAL_CHARS = 1000
    
    # Base for yinyang wuxing symbols (2 states * 5 elements)
    YIN_YANG_WUXING_BASE = 10
    # Base for bagua symbols (8 trigrams)
    BAGUA_BASE = 8

    # Symbols needed to encode ~1000 characters
    symbols_per_char_circle = math.ceil(math.log(TOTAL_CHARS, YIN_YANG_WUXING_BASE))
    symbols_per_char_square = math.ceil(math.log(TOTAL_CHARS, BAGUA_BASE))

    # Total characters that can be engraved on each plate
    chars_per_circle = math.floor(SYMBOLS_PER_CIRCLE / symbols_per_char_circle)
    chars_per_square = math.floor(SYMBOLS_PER_SQUARE / symbols_per_char_square)

    # --- Step 2 & 3: Determine the optimal packing of plates (Maximize circles) ---
    # We use a hexagonal packing of 3 rows (3-2-3) for the circles.
    circle_d = CIRCLE_R * 2
    
    # Width required for 3 circles in a row
    pack_width = 3 * circle_d
    # Vertical distance between the centers of hexagonally packed rows
    h_dist = math.sqrt(circle_d**2 - (circle_d / 2)**2)
    # Total height for 3 rows
    pack_height = circle_d + 2 * h_dist

    # The optimal number of circles (N) is 8, as this packing fits.
    N = 8
    
    # --- Step 4: Calculate the number of squares (M) in the remaining space ---
    # The 8-circle packing occupies a 120cm x 109.28cm bounding box.
    # We place it at one end of the 140x110cm material.

    # Region 1: The clear vertical strip left over on the side.
    rem_strip_w = RECT_W - pack_width
    rem_strip_h = RECT_H
    squares_from_strip = math.floor(rem_strip_w / SQUARE_S) * math.floor(rem_strip_h / SQUARE_S)

    # Region 2: The two side "pockets" next to the middle row of 2 circles.
    # These pockets are 20cm wide and 40cm high.
    pocket_w = circle_d / 2
    pocket_h = circle_d
    squares_from_pockets = 2 * (math.floor(pocket_w / SQUARE_S) * math.floor(pocket_h / SQUARE_S))

    # Total number of squared plates (M)
    M = squares_from_strip + squares_from_pockets

    # --- Step 5: Calculate the total number of characters (K) ---
    K = (N * chars_per_circle) + (M * chars_per_square)

    # --- Final Output ---
    print("To maximize the number of engraved characters, we determine the following:")
    print(f"\n1. Value of each plate:")
    print(f"   - A circular plate holds {chars_per_circle} characters.")
    print(f"   - A squared plate holds {chars_per_square} characters.")
    
    print(f"\n2. Optimal Plate Configuration:")
    print(f"   - Number of circular plates (N): {N}")
    print(f"   - Number of squared plates (M): {M}")

    print(f"\n3. Final Calculation for Total Characters (K):")
    print(f"   K = ({N} * {chars_per_circle}) + ({M} * {chars_per_square})")
    final_equation = f"   K = {N * chars_per_circle} + {M * chars_per_square} = {K}"
    print(final_equation)

    print("\n--------------------------------------------------")
    print("The optimal production is N:M:K")
    print(f"{N}:{M}:{K}")
    print("--------------------------------------------------")

solve_laozi_engraving()
<<<8:38:30084>>>