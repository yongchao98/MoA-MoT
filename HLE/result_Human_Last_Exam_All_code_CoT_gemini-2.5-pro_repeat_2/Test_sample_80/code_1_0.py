def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares on an 8x8 chessboard
    after placing the maximum number of non-attacking bishops.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An N x N board has 4 sides with N-1 unique squares each.
    total_edge_squares = 4 * (board_size - 1)

    # Step 2: Calculate the maximum number of non-attacking bishops on the edge.
    # For an even-sized N x N board, we can place N-1 bishops on the white
    # edge squares and N-1 on the black edge squares.
    # For N=8, this is 7 on white and 7 on black.
    max_bishops_on_edge = 2 * (board_size - 1)

    # Step 3: Calculate the number of edge squares that will lack bishops.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    print("--- Chess Bishop Puzzle ---")
    print(f"Total number of edge squares on a {board_size}x{board_size} board: {total_edge_squares}")
    print(f"Maximum non-attacking bishops that can be placed on the edge: {max_bishops_on_edge}")
    print("\nTo find the number of squares lacking bishops, we calculate:")
    print(f"{total_edge_squares} (total edge squares) - {max_bishops_on_edge} (placed bishops) = {empty_edge_squares}")
    print(f"\nFinal Answer: {empty_edge_squares} edge squares would lack bishops.")

solve_bishop_puzzle()