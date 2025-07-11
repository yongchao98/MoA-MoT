def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing
    the maximum number of non-attacking bishops on the edge.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An N x N board has 4 sides, but the 4 corners are shared.
    # The formula is 4 * N - 4, or 4 * (N - 1).
    total_edge_squares = 4 * (board_size - 1)
    print(f"A standard {board_size}x{board_size} chessboard has {total_edge_squares} edge squares.")

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # The 28 edge squares are divided into 14 light and 14 dark squares.
    # Bishops on light squares don't attack bishops on dark squares.
    # Max bishops on 14 light edge squares = 7
    # Max bishops on 14 dark edge squares = 7
    # Total max bishops = 7 + 7 = 14.
    max_bishops = 2 * (board_size - 1)
    print(f"The maximum number of non-attacking bishops that can be placed on the edge is {max_bishops}.")

    # Step 3: Calculate the number of edge squares without bishops.
    empty_edge_squares = total_edge_squares - max_bishops
    print("The number of edge squares that would lack bishops is the total number of edge squares minus the maximum number of bishops.")

    # Step 4: Display the final equation and result.
    print(f"Final calculation: {total_edge_squares} - {max_bishops} = {empty_edge_squares}")

# Run the solver
solve_bishops_puzzle()
<<<14>>>