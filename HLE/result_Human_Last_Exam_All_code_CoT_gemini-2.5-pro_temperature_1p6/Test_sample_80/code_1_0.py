def solve_bishops_problem():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_side = 8

    # Step 1: Calculate the total number of edge squares on an 8x8 board.
    # An n x n board has 4 sides of n squares, but the 4 corners are counted
    # twice, so the formula is 4 * n - 4.
    total_edge_squares = 4 * board_side - 4
    print(f"A standard {board_side}x{board_side} chessboard has {total_edge_squares} edge squares.")
    print(f"The calculation is: 4 * {board_side} - 4 = {total_edge_squares}")
    print("-" * 50)

    # Step 2: Determine the maximum number of non-attacking bishops.
    # We can analyze the white and black squares separately.
    # There are 28 / 2 = 14 white edge squares and 14 black edge squares.
    white_edge_squares = total_edge_squares // 2
    black_edge_squares = total_edge_squares // 2

    print(f"To place the maximum number of bishops, they cannot attack each other.")
    print(f"Since bishops only attack along diagonals of the same color, we can analyze the {white_edge_squares} white squares and {black_edge_squares} black squares independently.")
    
    # By analyzing the diagonal attack patterns of the 14 white edge squares,
    # it can be shown that a maximum of 7 bishops can be placed.
    # A valid placement is on squares b1, d1, f1, h1, c8, e8, and g8.
    max_bishops_on_white = 7
    print(f"The maximum number of non-attacking bishops on the {white_edge_squares} white edge squares is {max_bishops_on_white}.")

    # By symmetry, the same logic applies to the black squares.
    max_bishops_on_black = 7
    print(f"By symmetry, the maximum number on the {black_edge_squares} black edge squares is also {max_bishops_on_black}.")
    print("-" * 50)

    # Step 3: Calculate the total number of placed bishops.
    total_bishops_placed = max_bishops_on_white + max_bishops_on_black
    print("The total number of bishops successfully placed is the sum of the two:")
    print(f"Total bishops = {max_bishops_on_white} (white) + {max_bishops_on_black} (black) = {total_bishops_placed}")
    print("-" * 50)

    # Step 4: Calculate the number of empty edge squares.
    empty_edge_squares = total_edge_squares - total_bishops_placed
    print("The number of edge squares that lack bishops is the total number of edge squares minus the bishops we placed.")
    print(f"Final Equation: {total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}")
    
    return empty_edge_squares

# Run the function to display the steps and get the answer.
final_answer = solve_bishops_problem()