def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of an 8x8 chessboard.
    """
    board_size = 8

    # 1. Calculate the total number of edge squares.
    # An 8x8 board has 8 squares on each of its 4 sides.
    # The 4 corner squares are counted twice, so we subtract 4.
    # Total edge squares = (8 * 4) - 4 = 28.
    # Alternatively: 8 (top) + 8 (bottom) + 6 (left side inner) + 6 (right side inner) = 28
    total_edge_squares = (board_size * 4) - 4

    # 2. Determine the maximum number of non-attacking bishops on edge squares.
    # Bishops attack diagonally. The maximum number of non-attacking bishops
    # that can be placed on an 8x8 board is 14.
    # A valid configuration for the edge squares is to place bishops on all 8 squares of
    # the 'a' file and all 8 squares of the 'h' file. This is 16 squares.
    # However, in this setup, the bishop on a1 attacks h8, and the bishop on a8 attacks h1.
    # To resolve these two conflicts, we must remove two bishops (e.g., h1 and h8).
    # This leaves 16 - 2 = 14 bishops. This is the maximum possible.
    max_bishops_on_edge = 14

    # 3. Calculate the number of empty edge squares.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    print(f"A standard chessboard has {total_edge_squares} edge squares.")
    print(f"The maximum number of non-attacking bishops that can be placed on them is {max_bishops_on_edge}.")
    print("To find the number of empty squares, we subtract the number of bishops from the total number of edge squares.")
    print(f"Final Equation: {total_edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")

solve_chess_puzzle()
<<<14>>>