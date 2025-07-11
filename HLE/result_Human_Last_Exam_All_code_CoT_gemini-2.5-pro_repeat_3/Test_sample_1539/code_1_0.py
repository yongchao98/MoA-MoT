import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and square plates to maximize
    the number of engraved characters.
    """
    # --- Parameters ---
    # Material dimensions
    MATERIAL_W = 140
    MATERIAL_H = 110

    # Circular plates ("Book of Heaven")
    CIRCLE_RADIUS = 20
    CIRCLE_DIAMETER = 40
    CIRCLE_CHARS = 9999

    # Square plates ("Book of Earth")
    SQUARE_SIDE = 10
    SQUARE_CHARS = 360

    # To store the best result found: [N, M, K]
    best_result = [0, 0, 0]

    # --- Strategy 1: Squares Only ---
    m_w = MATERIAL_W // SQUARE_SIDE
    m_h = MATERIAL_H // SQUARE_SIDE
    m_total = m_w * m_h
    k_total = m_total * SQUARE_CHARS
    
    if k_total > best_result[2]:
        best_result = [0, m_total, k_total]

    # --- Strategy 2: Grid Pack Circles + Squares from leftovers ---
    # We test both orientations of the material (e.g., 140x110 and 110x140)
    for sheet_w, sheet_h in [(MATERIAL_W, MATERIAL_H), (MATERIAL_H, MATERIAL_W)]:
        # Max circles in a simple grid
        n_w = sheet_w // CIRCLE_DIAMETER
        n_h = sheet_h // CIRCLE_DIAMETER
        n_circles = n_w * n_h
        
        # Area used by circles
        used_w = n_w * CIRCLE_DIAMETER
        used_h = n_h * CIRCLE_DIAMETER
        
        # Calculate squares from leftover rectangular strips
        m_squares = 0
        # Strip 1:
        rem_w1 = sheet_w - used_w
        rem_h1 = sheet_h
        m_squares += (rem_w1 // SQUARE_SIDE) * (rem_h1 // SQUARE_SIDE)
        
        # Strip 2:
        rem_w2 = used_w
        rem_h2 = sheet_h - used_h
        m_squares += (rem_w2 // SQUARE_SIDE) * (rem_h2 // SQUARE_SIDE)
        
        k_total = (n_circles * CIRCLE_CHARS) + (m_squares * SQUARE_CHARS)
        
        if k_total > best_result[2]:
            best_result = [n_circles, m_squares, k_total]

    # --- Strategy 3: Hexagonal Pack Circles (densest packing) ---
    # Assume M=0 as leftovers are irregular and likely unusable for squares.
    for sheet_w, sheet_h in [(MATERIAL_W, MATERIAL_H), (MATERIAL_H, MATERIAL_W)]:
        r = CIRCLE_RADIUS
        d = CIRCLE_DIAMETER
        
        # Check if the sheet is large enough for at least one circle
        if sheet_w < d or sheet_h < d:
            continue
            
        # Number of circles in the first and second (staggered) rows
        n_row1 = sheet_w // d
        n_row2 = (sheet_w - r) // d
        
        # Number of rows we can fit
        # Height of k rows = (k-1)*r*sqrt(3) + d
        # (num_rows - 1) * r * sqrt(3) + d <= sheet_h
        num_rows = int((sheet_h - d) / (r * math.sqrt(3))) + 1

        # Total circles is the sum of circles in alternating row types
        n_circles = math.ceil(num_rows / 2) * n_row1 + math.floor(num_rows / 2) * n_row2
        
        k_total = n_circles * CIRCLE_CHARS
        
        if k_total > best_result[2]:
            # For this strategy, we assume no squares can be salvaged (M=0)
            best_result = [n_circles, 0, k_total]

    # --- Final Result ---
    N, M, K = best_result
    print(f"{N}:{M}:{K}")

solve_laozi_engraving()
<<<9:0:89991>>>