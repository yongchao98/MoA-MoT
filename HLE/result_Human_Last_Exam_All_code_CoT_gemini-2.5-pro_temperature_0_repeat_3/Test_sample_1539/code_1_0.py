import math

def solve_laozi_puzzle():
    """
    Solves the optimization problem for Laozi's books.
    """
    SHEET_W, SHEET_H = 140, 110
    CIRCLE_D = 40
    SQUARE_S = 10
    CIRCLE_CHARS = 9999
    SQUARE_CHARS = 360

    best_n = 0
    best_m = 0
    max_k = 0

    def grid_pack(width, height, side):
        """Calculates how many squares/circles (using bounding box) fit in a grid."""
        if width < side or height < side:
            return 0
        return math.floor(width / side) * math.floor(height / side)

    def hex_pack(width, height, diameter):
        """Calculates max circles using hexagonal packing in two orientations."""
        if width < diameter or height < diameter:
            return 0
        
        radius = diameter / 2
        row_dist = diameter * math.sqrt(3) / 2

        # Orientation 1: Rows parallel to width
        n_per_row1 = math.floor(width / diameter)
        n_per_row2 = math.floor((width - radius) / diameter)
        if row_dist == 0:
            num_rows1 = 0
        else:
            num_rows1 = 1 + math.floor((height - diameter) / row_dist)
        
        if num_rows1 <= 0:
            count1 = 0
        else:
            count1 = math.ceil(num_rows1 / 2) * n_per_row1 + math.floor(num_rows1 / 2) * n_per_row2

        # Orientation 2: Rows parallel to height (swap width and height)
        n_per_row1 = math.floor(height / diameter)
        n_per_row2 = math.floor((height - radius) / diameter)
        if row_dist == 0:
            num_rows2 = 0
        else:
            num_rows2 = 1 + math.floor((width - diameter) / row_dist)

        if num_rows2 <= 0:
            count2 = 0
        else:
            count2 = math.ceil(num_rows2 / 2) * n_per_row1 + math.floor(num_rows2 / 2) * n_per_row2
        
        return max(count1, count2)

    # Strategy 1: Partitioning the sheet
    # Cut parallel to the 110cm side
    for w_c in range(0, SHEET_W + 1, 10):
        w_s = SHEET_W - w_c
        
        n = max(grid_pack(w_c, SHEET_H, CIRCLE_D), hex_pack(w_c, SHEET_H, CIRCLE_D))
        m = grid_pack(w_s, SHEET_H, SQUARE_S)
        k = n * CIRCLE_CHARS + m * SQUARE_CHARS
        
        if k > max_k:
            max_k = k
            best_n = n
            best_m = m

    # Cut parallel to the 140cm side
    for h_c in range(0, SHEET_H + 1, 10):
        h_s = SHEET_H - h_c
        
        n = max(grid_pack(SHEET_W, h_c, CIRCLE_D), hex_pack(SHEET_W, h_c, CIRCLE_D))
        m = grid_pack(SHEET_W, h_s, SQUARE_S)
        k = n * CIRCLE_CHARS + m * SQUARE_CHARS
        
        if k > max_k:
            max_k = k
            best_n = n
            best_m = m
            
    # Strategy 2: Full sheet for circles
    # Assume leftover scrap is unusable for 10x10 squares
    n_full_hex = hex_pack(SHEET_W, SHEET_H, CIRCLE_D)
    m_full_hex = 0
    k_full_hex = n_full_hex * CIRCLE_CHARS + m_full_hex * SQUARE_CHARS

    # Compare strategies and determine the final answer
    if k_full_hex > max_k:
        final_n = n_full_hex
        final_m = m_full_hex
        final_k = k_full_hex
    else:
        final_n = best_n
        final_m = best_m
        final_k = max_k

    # Print the final result
    print(f"The optimal production is {final_n} circular plates and {final_m} squared plates.")
    print("The calculation for the maximum number of characters is:")
    print(f"{final_n} * {CIRCLE_CHARS} + {final_m} * {SQUARE_CHARS} = {final_k}")
    
    # Hidden final answer for the system
    # print(f"<<<{final_n}:{final_m}:{final_k}>>>")

solve_laozi_puzzle()