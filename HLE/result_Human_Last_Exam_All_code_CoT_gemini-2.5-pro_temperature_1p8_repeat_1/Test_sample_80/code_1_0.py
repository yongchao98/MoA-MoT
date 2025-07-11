def solve_bishops_on_edge_puzzle():
    """
    Calculates the number of empty edge squares after placing the
    maximum number of non-attacking bishops on the edge of a chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of squares on the edge of the board.
    # An n x n board has 4 sides of (n-1) squares each.
    total_edge_squares = (board_size - 1) * 4

    # Step 2: Calculate the maximum number of non-attacking bishops on the edge.
    # Bishops on light squares don't attack bishops on dark squares. The 28 edge
    # squares are split into 14 light and 14 dark squares.
    # Edge squares of a single color form diagonal pairs. For each pair,
    # we can only place one bishop. So, the number of bishops is half the
    # number of edge squares.
    max_bishops_placed = total_edge_squares // 2

    # Step 3: Calculate how many edge squares lack bishops.
    # This is the total number of edge squares minus the number occupied by bishops.
    empty_edge_squares = total_edge_squares - max_bishops_placed

    # Step 4: Print the reasoning and the final calculation.
    print(f"A chessboard has {board_size}x{board_size} squares.")
    print(f"The number of edge squares is: ({board_size} - 1) * 4 = {total_edge_squares}")
    print(f"To place the maximum number of non-attacking bishops, we can place one bishop for every two edge squares. This gives us {max_bishops_placed} bishops.")
    print("The number of empty edge squares is the total number of edge squares minus the number of bishops placed.")
    print(f"Final equation: {total_edge_squares} - {max_bishops_placed} = {empty_edge_squares}")
    print(f"Therefore, {empty_edge_squares} edge squares would lack bishops.")

solve_bishops_on_edge_puzzle()
<<<14>>>