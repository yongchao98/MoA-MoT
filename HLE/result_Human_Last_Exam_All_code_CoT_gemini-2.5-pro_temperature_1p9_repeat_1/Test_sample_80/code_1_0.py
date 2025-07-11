def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on them.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    total_edge_squares = (board_size * 4) - 4

    print(f"A standard {board_size}x{board_size} chessboard has {total_edge_squares} edge squares.")

    # Step 2: Divide the problem by color.
    # The edge squares are split evenly between white and black.
    squares_per_color = total_edge_squares // 2
    print(f"These squares are divided into {squares_per_color} white and {squares_per_color} black squares.")

    # Step 3: Calculate the maximum bishops that can be placed.
    # On the edge, every square of a given color forms an attacking pair with
    # exactly one other edge square of the same color.
    # To place the maximum number of non-attacking bishops, we place one per pair.
    max_bishops_per_color = squares_per_color // 2

    # The total number of bishops is the sum for both colors.
    total_placed_bishops = max_bishops_per_color * 2

    print(f"To place the maximum number of non-attacking bishops, you can place {total_placed_bishops} bishops on the edge.")

    # Step 4: Calculate the number of empty squares.
    empty_squares = total_edge_squares - total_placed_bishops

    print("\nThe number of edge squares that would lack bishops is calculated as follows:")
    print(f"{total_edge_squares} (total edge squares) - {total_placed_bishops} (placed bishops) = {empty_squares} (empty squares)")

solve_bishop_puzzle()
<<<14>>>