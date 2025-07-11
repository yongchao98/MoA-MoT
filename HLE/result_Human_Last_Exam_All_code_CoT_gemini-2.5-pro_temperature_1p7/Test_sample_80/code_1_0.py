def solve_bishops_puzzle():
    """
    Calculates the number of edge squares lacking bishops when the maximum number
    of non-attacking bishops are placed on the edge of an 8x8 chessboard.
    """
    # Step 1: Define board size and calculate total edge squares.
    board_size = 8
    total_squares = board_size * board_size
    inner_squares = (board_size - 2) * (board_size - 2)
    edge_squares = total_squares - inner_squares

    print(f"A standard {board_size}x{board_size} chessboard has {edge_squares} edge squares.")
    print("These edge squares are composed of 14 light-colored squares and 14 dark-colored squares.")
    print("-" * 20)

    # Step 2: Determine the maximum number of bishops that can be placed.
    # Bishops on light squares don't attack bishops on dark squares.
    # The 14 squares of each color form 7 pairs that lie on the same diagonals.
    # Therefore, we can only place one bishop for each of these 7 diagonals for each color.
    max_bishops_per_color = 7
    total_bishops = max_bishops_per_color * 2

    print(f"On the 14 light squares, a maximum of {max_bishops_per_color} non-attacking bishops can be placed.")
    print(f"On the 14 dark squares, a maximum of {max_bishops_per_color} non-attacking bishops can be placed.")
    print(f"The total number of bishops placed on the edge is {max_bishops_per_color} + {max_bishops_per_color} = {total_bishops}.")
    print("-" * 20)

    # Step 3: Calculate the number of empty squares.
    lacking_bishops = edge_squares - total_bishops

    print("The number of edge squares that would lack bishops is calculated as follows:")
    # The final print statement shows the equation with the actual numbers.
    print(f"{edge_squares} (total edge squares) - {total_bishops} (placed bishops) = {lacking_bishops}")

solve_bishops_puzzle()
<<<14>>>