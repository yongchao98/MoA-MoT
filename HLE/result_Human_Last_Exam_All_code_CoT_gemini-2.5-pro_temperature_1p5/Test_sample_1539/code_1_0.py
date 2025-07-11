import math

def solve_laozi_engraving():
    """
    Solves the optimization problem for Laozi's books by testing different packing strategies.
    """
    # --- Problem Constants ---
    SHEET_W, SHEET_H = 140.0, 110.0
    CIRC_D, CIRC_R = 40.0, 20.0
    SQ_S = 10.0
    CHAR_PER_CIRC = 9999
    CHAR_PER_SQ = 360

    best_N, best_M, best_K = 0, 0, 0
    best_strategy = "None"

    # --- Strategy 1: Squares Only ---
    # Use the entire 140x110 sheet for 10x10 squares.
    n1 = 0
    m1 = math.floor(SHEET_W / SQ_S) * math.floor(SHEET_H / SQ_S)
    k1 = n1 * CHAR_PER_CIRC + m1 * CHAR_PER_SQ
    if k1 > best_K:
        best_N, best_M, best_K, best_strategy = n1, m1, k1, "Squares Only"

    # --- Strategy 2: Grid-Packed Circles + Squares ---
    # Pack circles in a simple grid and use the leftover L-shaped area for squares.
    n2_cols = math.floor(SHEET_W / CIRC_D)  # 140/40 = 3
    n2_rows = math.floor(SHEET_H / CIRC_D)  # 110/40 = 2
    n2 = n2_cols * n2_rows # 3 * 2 = 6 circles

    used_w = n2_cols * CIRC_D # 120cm
    used_h = n2_rows * CIRC_D # 80cm
    
    # Leftover consists of two rectangular strips
    rem_strip1_w, rem_strip1_h = SHEET_W, SHEET_H - used_h      # 140 x 30
    rem_strip2_w, rem_strip2_h = SHEET_W - used_w, used_h      # 20 x 80
    
    m2_1 = math.floor(rem_strip1_w / SQ_S) * math.floor(rem_strip1_h / SQ_S) # 14*3=42
    m2_2 = math.floor(rem_strip2_w / SQ_S) * math.floor(rem_strip2_h / SQ_S) # 2*8=16
    m2 = m2_1 + m2_2 # 42 + 16 = 58 squares

    k2 = n2 * CHAR_PER_CIRC + m2 * CHAR_PER_SQ
    if k2 > best_K:
        best_N, best_M, best_K, best_strategy = n2, m2, k2, "Grid-Packed Circles"

    # --- Strategy 3: Hexagonally-Packed Circles + Squares in Scraps ---
    # This is the most efficient packing for circles.
    # Height of a staggered row: r * sqrt(3)
    h_stagger = CIRC_R * math.sqrt(3) # ~34.64 cm
    
    # We found that 3 rows of 3 circles fit in a 140x110 sheet (height: 40 + 2*h_stagger = 109.28 cm)
    n3 = 9 
    
    # This packing leaves small rectangular areas on the sides.
    # We calculate the number of squares that can fit in these scraps.
    # Scrap A (bottom-right): 20cm width, h_stagger height
    m3_a = math.floor(20.0 / SQ_S) * math.floor(h_stagger / SQ_S) # 2 * 3 = 6
    # Scrap B (top-right): 20cm width, 110 - (top edge of 2nd row) height
    y_r2_max = CIRC_D + h_stagger
    m3_b = math.floor(20.0 / SQ_S) * math.floor((SHEET_H - y_r2_max) / SQ_S) # 2 * 3 = 6
    # Scrap C (middle-left): 20cm width, (bottom edge of 3rd row - top edge of 1st row) height
    y_r3_min = CIRC_D + h_stagger - CIRC_R
    y_r1_max = CIRC_D
    m3_c = math.floor(20.0 / SQ_S) * math.floor((y_r3_min - y_r1_max) / SQ_S) # 2 * 2 = 4
    
    m3 = m3_a + m3_b + m3_c # 6 + 6 + 4 = 16 squares
    
    k3 = n3 * CHAR_PER_CIRC + m3 * CHAR_PER_SQ
    if k3 > best_K:
        best_N, best_M, best_K, best_strategy = n3, m3, k3, "Dense Hexagonal Packing"
        
    # --- Final Output ---
    print(f"The optimal strategy is '{best_strategy}', which maximizes the number of engraved characters.")
    print("\nThis method produces:")
    print(f"N = {best_N} circular plates")
    print(f"M = {best_M} squared plates")
    print("\nThe final calculation for the total characters (K) is:")
    print(f"{best_N} * {CHAR_PER_CIRC} + {best_M} * {CHAR_PER_SQ} = {best_N * CHAR_PER_CIRC + best_M * CHAR_PER_SQ}")
    print("\nAnswer in the required format N:M:K")
    print(f"{best_N}:{best_M}:{best_K}")

solve_laozi_engraving()