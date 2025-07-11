import math

def solve_laozi_puzzle():
    """
    Calculates the maximum number of Chinese characters that can be engraved
    by optimizing the cutting of circular and squared plates from a meteorite sheet.
    """

    # --- Step 1: Define Constants and Calculate Character Capacity ---

    # Material and plate dimensions in cm
    sheet_width = 140
    sheet_height = 110
    circle_radius = 20
    square_side = 10

    # Engraving capacity in symbols
    symbols_per_circle = 9999
    symbols_per_square = 360

    # Encoding calculation to find symbols per character
    # Book of Heaven (Yinyang Wuxing): 2 states (yin/yang) * 5 elements (wuxing) = 10 unique symbols.
    # To represent ~1000 characters, we need 10^3 = 1000 combinations.
    symbols_per_char_heaven = 3
    # Book of Earth (Bagua): 8 unique trigrams.
    # To represent ~1000 characters, 8^3=512 is not enough, so we need 8^4 = 4096 combinations.
    symbols_per_char_earth = 4

    # Calculate characters per plate
    chars_per_circle = symbols_per_circle // symbols_per_char_heaven
    chars_per_square = symbols_per_square // symbols_per_char_earth

    print("Step-by-step thinking process:")
    print("1. Determine the number of symbols needed to encode one Chinese character.")
    print(f"   - For the Book of Heaven, there are 10 unique yinyang wuxing symbols. To represent ~1000 characters, {symbols_per_char_heaven} symbols are needed per character (10^{symbols_per_char_heaven} >= 1000).")
    print(f"   - Characters per circular plate = {symbols_per_circle} symbols / {symbols_per_char_heaven} symbols/char = {chars_per_circle} characters.")
    print(f"   - For the Book of Earth, there are 8 unique bagua symbols. To represent ~1000 characters, {symbols_per_char_earth} symbols are needed per character (8^{symbols_per_char_earth-1} < 1000, 8^{symbols_per_char_earth} >= 1000).")
    print(f"   - Characters per squared plate = {symbols_per_square} symbols / {symbols_per_char_earth} symbols/char = {chars_per_square} characters.")

    # --- Step 2: Determine Optimal Packing Strategy ---

    # We must maximize circles first. Using a hexagonal packing is more efficient than a simple grid.
    # A 3-2-3 arrangement of 8 circles (radius 20cm) fits within a 120cm x 109.28cm area.
    # This fits inside the 140cm x 110cm sheet.
    num_circles = 8
    area_for_circles_w = 120
    
    print("\n2. Determine the optimal number of plates to cut.")
    print("   - To maximize characters, we must prioritize cutting the circular plates.")
    print(f"   - Using an efficient hexagonal packing, we can fit {num_circles} circular plates.")
    print(f"   - This arrangement occupies a {area_for_circles_w}cm x ~109.3cm section of the {sheet_width}cm x {sheet_height}cm material.")
    print(f"   - The remaining material is a strip of ({sheet_width} - {area_for_circles_w})cm x {sheet_height}cm = 20cm x 110cm.")
    
    # Calculate how many squares fit in the remaining rectangular strip.
    rem_width = sheet_width - area_for_circles_w
    rem_height = sheet_height
    num_squares_in_rem_w = rem_width // square_side
    num_squares_in_rem_h = rem_height // square_side
    num_squares = num_squares_in_rem_w * num_squares_in_rem_h
    
    print(f"   - Number of squared plates (M) from this strip = floor({rem_width}/{square_side}) * floor({rem_height}/{square_side}) = {num_squares_in_rem_w} * {num_squares_in_rem_h} = {num_squares}.")

    # --- Step 3: Calculate Total Characters and Format Output ---
    N = num_circles
    M = num_squares
    K = (N * chars_per_circle) + (M * chars_per_square)

    print("\n3. Calculate the maximum number of engraved characters (K).")
    print("   The final equation is: K = N * (chars per circle) + M * (chars per square).")
    print(f"   K = {N} * {chars_per_circle} + {M} * {chars_per_square} = {N * chars_per_circle} + {M * chars_per_square} = {K}")

    final_answer = f"{N}:{M}:{K}"
    print(f"\nThe optimal result is N={N} circular plates, M={M} squared plates, and K={K} total characters.")
    
    print(f"<<<{final_answer}>>>")

solve_laozi_puzzle()