def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    total_squares = board_size * board_size
    inner_squares = (board_size - 2) * (board_size - 2)
    edge_squares = total_squares - inner_squares

    print(f"An {board_size}x{board_size} chessboard has {edge_squares} edge squares.")

    # Step 2: Determine the number of squares per color on the edge.
    # On an empty board, the edge is composed of an equal number of light and dark squares.
    edge_squares_per_color = edge_squares // 2

    print(f"These {edge_squares} edge squares are made up of {edge_squares_per_color} light squares and {edge_squares_per_color} dark squares.")

    # Step 3: Determine the maximum number of non-attacking bishops.
    # Bishops only move on diagonals of their starting color.
    # To place the maximum number of non-attacking bishops, we can occupy all edge squares of a single color.
    max_bishops_on_edge = edge_squares_per_color

    print(f"The maximum number of non-attacking bishops that can be placed on the edge is {max_bishops_on_edge}.")

    # Step 4: Calculate the number of empty edge squares.
    empty_edge_squares = edge_squares - max_bishops_on_edge

    print("\nTo find the number of edge squares that lack bishops, we subtract the number of placed bishops from the total number of edge squares.")
    print(f"The final calculation is: {edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")


solve_chess_puzzle()
<<<14>>>