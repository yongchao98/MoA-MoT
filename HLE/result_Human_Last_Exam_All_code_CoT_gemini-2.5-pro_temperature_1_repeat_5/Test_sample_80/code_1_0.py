def solve_bishops_on_edge_problem():
    """
    Calculates the number of empty edge squares on an 8x8 chessboard
    after placing the maximum number of non-attacking bishops.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An n x n board has 4*(n-1) squares on its edge.
    total_edge_squares = 4 * (board_size - 1)

    # Step 2: Determine the maximum number of non-attacking bishops
    # that can be placed on the edge of an n x n board.
    # This is 2 for each color for the corners and n-2 for each side,
    # leading to the formula 2*n - 2.
    # For an 8x8 board, this is 2 * 8 - 2 = 14.
    max_bishops_on_edge = 2 * board_size - 2

    # Step 3: Calculate the number of edge squares without bishops.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    # Print the explanation and the final equation with the calculated numbers.
    print(f"A standard {board_size}x{board_size} chessboard has {total_edge_squares} edge squares.")
    print(f"The maximum number of non-attacking bishops that can be placed on these squares is {max_bishops_on_edge}.")
    print("Therefore, the number of edge squares that would lack bishops is:")
    print(f"{total_edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")

solve_bishops_on_edge_problem()
<<<14>>>