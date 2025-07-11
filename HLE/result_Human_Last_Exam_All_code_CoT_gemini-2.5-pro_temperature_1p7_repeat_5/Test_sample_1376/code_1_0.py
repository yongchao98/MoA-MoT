import math

def solve_tomb_engraving():
    """
    Calculates the optimal number of squares (N) and circles (M)
    to maximize the total number of engraved characters (K).
    """
    # --- Step 1: Define constants and problem parameters ---

    # Material dimensions in cm
    sheet_w, sheet_h = 140, 110

    # Artifact dimensions in cm
    square_side = 10
    # A circle with 20cm radius requires a 40x40cm bounding box
    circle_bounding_box_side = 40

    # Characters per artifact
    chars_per_square = 4
    # Calculate characters per circle based on Bagua encoding
    # Bio uses 7x7x9 = 441 unique characters
    # Bagua has 8 symbols. We need log8(441) symbols per character.
    # math.ceil(math.log(441, 8)) = 3 symbols per character.
    symbols_per_character = 3
    symbols_per_circle = 999
    chars_per_circle = math.floor(symbols_per_circle / symbols_per_character)

    # --- Step 2: Maximize the number of circles (M) ---
    # The strategy is to fit as many circles as possible, as they hold more characters.
    # We use integer division to find how many 40x40 boxes fit.

    # Option 1: Sheet as 140x110
    m_cols_1 = sheet_w // circle_bounding_box_side
    m_rows_1 = sheet_h // circle_bounding_box_side
    m_option_1 = m_cols_1 * m_rows_1

    # Option 2: Sheet as 110x140 (rotated)
    m_cols_2 = sheet_h // circle_bounding_box_side
    m_rows_2 = sheet_w // circle_bounding_box_side
    m_option_2 = m_cols_2 * m_rows_2
    
    # We choose the orientation that fits more circles (in this case, they are equal).
    if m_option_1 >= m_option_2:
        M = m_option_1
        m_cols = m_cols_1
        m_rows = m_rows_1
    else:
        M = m_option_2
        m_cols = m_cols_2
        m_rows = m_rows_2
        # Swap sheet dimensions to match the chosen orientation for calculation
        sheet_w, sheet_h = sheet_h, sheet_w
    
    # --- Step 3: Calculate remaining area and the number of squares (N) ---
    
    # Area used by the circles' bounding boxes
    used_width = m_cols * circle_bounding_box_side
    used_height = m_rows * circle_bounding_box_side

    # The remaining L-shaped area is split into two rectangles
    # Rectangle 1
    rem_rect1_w = sheet_w - used_width
    rem_rect1_h = sheet_h
    n1 = (rem_rect1_w // square_side) * (rem_rect1_h // square_side)

    # Rectangle 2
    rem_rect2_w = used_width
    rem_rect2_h = sheet_h - used_height
    n2 = (rem_rect2_w // square_side) * (rem_rect2_h // square_side)

    N = n1 + n2
    
    # --- Step 4: Calculate the total number of characters (K) ---
    K = (N * chars_per_square) + (M * chars_per_circle)

    # --- Step 5: Print the results ---
    print("--- Maximization Plan ---")
    print(f"Optimal number of circles (M) to maximize characters: {M}")
    print(f"Optimal number of squares (N) from remaining material: {N}")
    print("\n--- Final Calculation ---")
    print("K = (N * characters_per_square) + (M * characters_per_circle)")
    print(f"K = ({N} * {chars_per_square}) + ({M} * {chars_per_circle})")
    n_chars = N * chars_per_square
    m_chars = M * chars_per_circle
    print(f"K = {n_chars} + {m_chars}")
    print(f"K = {K}")
    
    print("\n--- Answer ---")
    print(f"The final answer in N:M:K format is {N}:{M}:{K}")


# Execute the function to find the solution
solve_tomb_engraving()