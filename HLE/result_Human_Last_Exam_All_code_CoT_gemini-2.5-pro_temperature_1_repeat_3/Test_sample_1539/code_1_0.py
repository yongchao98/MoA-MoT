import math

def solve():
    """
    Calculates the maximum number of Chinese characters that can be engraved
    by optimizing the cutting of circular and squared plates from a meteorite sheet.
    """

    # --- Step 1: Define constants ---
    sheet_width = 140  # cm
    sheet_height = 110 # cm
    
    # Circular plate ("Book of Heaven")
    circle_radius = 20
    circle_diameter = 2 * circle_radius
    circle_symbols = 9999
    
    # Squared plate ("Book of Earth")
    square_side = 10
    square_symbols = 360
    
    # Encoding parameters
    num_chinese_chars = 1000
    yinyang_states = 2
    wuxing_states = 5
    bagua_states = 8

    # --- Step 2: Calculate encoding efficiency and character value per plate ---
    
    # For circular plates (yinyang wuxing)
    yinyang_wuxing_combinations = yinyang_states * wuxing_states
    symbols_per_char_circle = math.ceil(math.log(num_chinese_chars, yinyang_wuxing_combinations))
    chars_per_circle = circle_symbols // symbols_per_char_circle
    
    # For squared plates (bagua)
    symbols_per_char_square = math.ceil(math.log(num_chinese_chars, bagua_states))
    chars_per_square = square_symbols // symbols_per_char_square

    # --- Step 3: Solve the packing problem ---
    # The value of a circle (3333 chars) is much higher than a square (90 chars).
    # The optimal strategy is to maximize the number of circles.
    
    # We analyze packing 40cm-diameter circles in a 140x110cm rectangle.
    # A staggered (hexagonal) layout is most efficient.
    # Let's align the rows with the 140cm side.
    
    # A row can fit floor(140/40) = 3 circles.
    # A staggered row can also fit 3 circles, occupying the full 140cm width.
    # The vertical distance between centers of staggered rows is D * sqrt(3)/2.
    row_center_dist_v = circle_diameter * math.sqrt(3) / 2
    
    # Let's find how many rows (k) can fit in the 110cm height.
    # Total height for k rows = (k-1) * row_center_dist_v + circle_diameter
    # (k-1) * 34.64 + 40 <= 110  => (k-1) <= 70 / 34.64 => k-1 <= 2.02 => k <= 3.02
    num_rows = 3
    
    # With 3 staggered rows of 3 circles each, we can fit 3 * 3 = 9 circles.
    # This layout (3-3-3 staggered) occupies 140cm in width and
    # (3-1)*34.64 + 40 = 109.28cm in height. This fits within the 140x110cm sheet.
    num_circles = 9
    
    # This packing leaves no significant material to cut 10x10cm squares.
    # The remaining area is a thin 140cm x 0.72cm strip and small interstitial
    # spaces between circles, from which no squares can be cut.
    num_squares = 0
    
    # --- Step 4: Calculate the total number of characters ---
    max_chars = (num_circles * chars_per_circle) + (num_squares * chars_per_square)

    # --- Step 5: Print the final output ---
    print("This problem is an optimization puzzle. The goal is to maximize the total number of engraved characters.")
    print("First, we calculate the number of characters per plate:")
    print(f"  - A circular plate encodes {chars_per_circle} characters ({circle_symbols} symbols / {symbols_per_char_circle} symbols per character).")
    print(f"  - A squared plate encodes {chars_per_square} characters ({square_symbols} symbols / {symbols_per_char_square} symbols per character).")
    print("\nSince circular plates are far more valuable, we prioritize cutting the maximum number of circles.")
    print("Using an efficient staggered packing, we can cut 9 circles from the 140x110cm sheet.")
    print("This leaves no usable material for squared plates.\n")
    
    print("The final calculation is:")
    # The final equation as requested.
    print(f"{num_circles} circular plates * {chars_per_circle} characters/plate + {num_squares} squared plates * {chars_per_square} characters/plate = {max_chars} characters")
    
    final_answer = f"{num_circles}:{num_squares}:{max_chars}"
    print(f"\n<<<N:M:K = {final_answer}>>>")


solve()
<<<9:0:29997>>>