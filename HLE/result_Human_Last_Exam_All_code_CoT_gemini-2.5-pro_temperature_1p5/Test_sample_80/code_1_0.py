def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing
    the maximum possible number of non-attacking bishops on them.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares on an N x N board.
    # Formula: 4 * (N - 1)
    total_edge_squares = 4 * (board_size - 1)

    # Step 2: Determine the maximum number of non-attacking bishops.
    # An 8x8 board has 28 edge squares (14 white, 14 black).
    # Bishops on squares of different colors never attack each other.
    # For each color, the 14 edge squares form 7 pairs of squares that
    # lie on the same diagonal. From each pair, only one bishop can be placed.
    # Therefore, we can place 7 bishops on white squares and 7 on black squares.
    max_bishops_on_white_edges = 7
    max_bishops_on_black_edges = 7
    total_placed_bishops = max_bishops_on_white_edges + max_bishops_on_black_edges

    # Step 3: Calculate how many edge squares are empty.
    empty_edge_squares = total_edge_squares - total_placed_bishops

    # Step 4: Print the final equation and the result.
    print(f"A standard {board_size}x{board_size} chessboard has {total_edge_squares} edge squares.")
    print(f"The maximum number of non-attacking bishops that can be placed on these squares is {total_placed_bishops}.")
    print("The number of empty edge squares is the difference between these two numbers.")
    print("Final Equation:")
    print(f"{total_edge_squares} - {total_placed_bishops} = {empty_edge_squares}")


solve_chess_puzzle()
<<<14>>>