import math

def solve_laozi_engraving():
    """
    Calculates the optimal number of circular and squared plates to maximize
    the number of engraved characters from a 140x110cm sheet.
    """
    
    # --- Constants ---
    SHEET_W = 140
    SHEET_H = 110
    CIRCLE_DIAMETER = 40
    SQUARE_SIDE = 10
    CIRCLE_CHARS = 9999
    SQUARE_CHARS = 360

    # --- Scenario 1: Only Squares ---
    n1 = 0
    m1 = math.floor(SHEET_W / SQUARE_SIDE) * math.floor(SHEET_H / SQUARE_SIDE)
    k1 = n1 * CIRCLE_CHARS + m1 * SQUARE_CHARS

    # --- Scenario 2: Grid Packing of 6 Circles ---
    # A 3x2 grid of 40cm circles uses 120x80cm.
    n2 = 6
    # Remaining L-shape can be split into a 140x30cm rectangle and a 20x80cm rectangle.
    m2_rect1 = math.floor(140 / SQUARE_SIDE) * math.floor(30 / SQUARE_SIDE)
    m2_rect2 = math.floor(20 / SQUARE_SIDE) * math.floor(80 / SQUARE_SIDE)
    m2 = m2_rect1 + m2_rect2
    k2 = n2 * CIRCLE_CHARS + m2 * SQUARE_CHARS

    # --- Scenario 3: Staggered Packing of 8 Circles ---
    # A 3-2-3 staggered packing fits in a ~120x110cm bounding box.
    n3 = 8
    # This leaves a main rectangular strip of (140-120) x 110 cm.
    m3 = math.floor(20 / SQUARE_SIDE) * math.floor(110 / SQUARE_SIDE)
    k3 = n3 * CIRCLE_CHARS + m3 * SQUARE_CHARS

    # --- Scenario 4: Maximum Staggered Packing of 9 Circles ---
    # A 3x3 staggered packing of 9 circles fits in a 140cm x 109.28cm area.
    n4 = 9
    # The remaining area (140cm x 0.72cm) is too small for any 10x10cm squares.
    m4 = 0
    k4 = n4 * CIRCLE_CHARS + m4 * SQUARE_CHARS

    # --- Find the best scenario ---
    scenarios = {
        k1: (n1, m1),
        k2: (n2, m2),
        k3: (n3, m3),
        k4: (n4, m4)
    }

    best_k = max(scenarios.keys())
    best_n, best_m = scenarios[best_k]

    # --- Print the final answer and equation ---
    print("The optimal solution is to produce {} circular plates and {} squared plates.".format(best_n, best_m))
    print("The final equation for the maximum number of characters is:")
    print("{} * {} + {} * {} = {}".format(best_n, CIRCLE_CHARS, best_m, SQUARE_CHARS, best_k))
    print("\nFinal Answer in N:M:K format:")
    print("{}:{}:{}".format(best_n, best_m, best_k))


solve_laozi_engraving()
<<<9:0:89991>>>