def solve_alien_chess_colony():
    """
    Calculates the maximal size (K) of the aliens' colony.

    The solution is based on the property that the final shape of the colony
    is the smallest axis-aligned rectangle (bounding box) that contains all
    the initial 8 squares. The goal is to place the 6 non-fixed squares
    to maximize the area of this rectangle.
    """

    # 1. Define the board and coordinate system.
    # An 8x8 board corresponds to coordinates (0,0) to (7,7).
    # Let's map chess notation where a1 is (0,0), h1 is (0,7), a8 is (7,0), h8 is (7,7).
    # This means rows '1'-'8' map to 0-7 and columns 'a'-'h' map to 0-7.
    # The fixed squares are d5 and e5.
    # d5 -> row '5' is 4, col 'd' is 3. Coordinate: (4, 3)
    # e5 -> row '5' is 4, col 'e' is 4. Coordinate: (4, 4)
    fixed_squares = [(4, 3), (4, 4)]

    # 2. Define the boundaries of the final rectangle.
    # The rectangle is defined by r_min, r_max, c_min, c_max, which are the
    # minimum and maximum row and column indices of the initial 8 squares.
    # The fixed squares impose initial constraints on these boundaries.
    # r_min must be <= 4, and r_max must be >= 4.
    # c_min must be <= 3, and c_max must be >= 4.

    # 3. Optimize the placement of the 6 other squares.
    # To maximize the area of the rectangle, (r_max - r_min + 1) * (c_max - c_min + 1),
    # we need to make the bounding box as large as possible. The board itself
    # provides the ultimate limits: rows and columns are in the range [0, 7].

    # We can use our 6 disposable squares to push the boundaries to these limits.
    # - To make r_min as small as possible, we place a square in row 0.
    # - To make r_max as large as possible, we place a square in row 7.
    # - To make c_min as small as possible, we place a square in column 0.
    # - To make c_max as large as possible, we place a square in column 7.
    # For example, placing squares at the corners (0,0), (0,7), (7,0), and (7,7)
    # uses 4 of our 6 squares and guarantees the bounding box is the entire board.

    # 4. Calculate the dimensions and area of the maximal rectangle.
    r_min_optimal = 0
    r_max_optimal = 7
    c_min_optimal = 0
    c_max_optimal = 7

    height = r_max_optimal - r_min_optimal + 1
    width = c_max_optimal - c_min_optimal + 1

    K = height * width

    print("The maximal size K is the area of the largest possible bounding rectangle.")
    print(f"The optimal boundaries for the rectangle are:")
    print(f"  min row: {r_min_optimal}, max row: {r_max_optimal}")
    print(f"  min col: {c_min_optimal}, max col: {c_max_optimal}")
    print("\nCalculating the dimensions:")
    print(f"Height = {r_max_optimal} - {r_min_optimal} + 1 = {height}")
    print(f"Width  = {c_max_optimal} - {c_min_optimal} + 1 = {width}")
    print("\nThe final equation for the maximal area K is:")
    print(f"K = {height} * {width} = {K}")

solve_alien_chess_colony()