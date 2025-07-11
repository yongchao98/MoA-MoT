import math

def solve_laozi_puzzle():
    """
    Solves the optimization problem for Laozi's books.
    
    This function calculates the optimal number of circular and squared plates
    to maximize the total number of engraved characters, based on a specific
    packing strategy.
    """

    # --- Step 1: Define constants and calculate character capacity per plate ---
    
    # Material dimensions
    MATERIAL_W = 140  # cm
    MATERIAL_H = 110  # cm

    # Plate properties
    CIRCLE_RADIUS = 20
    CIRCLE_DIAMETER = CIRCLE_RADIUS * 2
    SQUARE_SIDE = 10

    # Symbol capacity
    SYMBOLS_PER_CIRCLE = 9999
    SYMBOLS_PER_SQUARE = 360
    
    # Encoding parameters
    NUM_CHARACTERS = 1000
    YINYANG_WUXING_ALPHABET_SIZE = 10 # (yin/yang * 5 elements)
    BAGUA_ALPHABET_SIZE = 8

    # Calculate characters per plate based on encoding principles
    # Symbols needed = ceil(log_alphabet(num_items))
    symbols_per_char_heaven = math.ceil(math.log(NUM_CHARACTERS, YINYANG_WUXING_ALPHABET_SIZE))
    chars_per_circle = math.floor(SYMBOLS_PER_CIRCLE / symbols_per_char_heaven)

    symbols_per_char_earth = math.ceil(math.log(NUM_CHARACTERS, BAGUA_ALPHABET_SIZE))
    chars_per_square = math.floor(SYMBOLS_PER_SQUARE / symbols_per_char_earth)

    print("--- Analysis ---")
    print(f"Book of Heaven (Circle): Requires {symbols_per_char_heaven} yinyang wuxing symbols per character. Each plate holds {chars_per_circle} characters.")
    print(f"Book of Earth (Square): Requires {symbols_per_char_earth} bagua symbols per character. Each plate holds {chars_per_square} characters.")
    print("-" * 16)

    # --- Step 2 & 3: Devise a packing strategy and calculate N and M ---
    
    # Strategy: Prioritize circles with hexagonal packing due to higher value, then fill leftovers with squares.
    # We can fit 3 rows of circles hexagonally within the 110cm height.
    # Vertical distance between rows of circles = sqrt(D^2 - (D/2)^2)
    row_dist_y = math.sqrt(CIRCLE_DIAMETER**2 - (CIRCLE_DIAMETER/2)**2)
    # Total height for 3 rows = Diameter + 2 * row_dist_y
    pack_height = CIRCLE_DIAMETER + 2 * row_dist_y # = 40 + 2*34.64 = 109.28 cm
    
    # Row arrangement for maximum circles:
    # Row 1 fits floor(140/40) = 3 circles.
    # Row 2 (offset) fits 2 circles.
    # Row 3 fits 3 circles.
    # This packing fits within a 120cm x 109.28cm bounding box.
    
    N_circles = 3 + 2 + 3
    pack_width = 3 * CIRCLE_DIAMETER - CIRCLE_DIAMETER # = 120 cm

    # The packing fits in the 140x110 material.
    # After packing the circles, a large rectangular piece of material is left over.
    leftover_w = MATERIAL_W - pack_width
    leftover_h = MATERIAL_H
    
    # Calculate how many squares fit in this leftover rectangular area.
    M_squares_in_leftover = math.floor(leftover_w / SQUARE_SIDE) * math.floor(leftover_h / SQUARE_SIDE)
    # The other small leftover areas (e.g., curved gaps between circles) cannot fit a 10x10cm square.

    print("--- Packing Strategy ---")
    print(f"Using a hexagonal layout, we can pack a maximum of {N_circles} circular plates.")
    print(f"This leaves a rectangular area of {leftover_w}cm x {leftover_h}cm.")
    print(f"From this leftover area, we can cut {M_squares_in_leftover} squared plates.")
    print("-" * 16)

    # --- Step 4: Calculate the maximal number of characters K ---
    K_max_chars = (N_circles * chars_per_circle) + (M_squares_in_leftover * chars_per_square)
    
    # --- Step 5: Output the result ---
    print("--- Final Calculation ---")
    final_equation = (
        f"{N_circles} circular plates * {chars_per_circle} chars/plate + "
        f"{M_squares_in_leftover} squared plates * {chars_per_square} chars/plate"
    )
    print(final_equation)
    print(f"= {N_circles * chars_per_circle} + {M_squares_in_leftover * chars_per_square}")
    print(f"= {int(K_max_chars)} characters")
    
    final_answer = f"{N_circles}:{M_squares_in_leftover}:{int(K_max_chars)}"
    print("\nFinal Answer (N:M:K):")
    print(final_answer)
    return final_answer

if __name__ == '__main__':
    solve_laozi_puzzle()