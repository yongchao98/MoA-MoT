def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_side = 8

    # Step 1: Calculate the total number of edge squares.
    # An n x n board has 4 sides of n squares, but the 4 corners are counted twice.
    # So, the formula is (n * 4) - 4.
    total_edge_squares = (board_side * 4) - 4
    print(f"A standard {board_side}x{board_side} board has {total_edge_squares} edge squares.")

    # Step 2: Determine the maximum number of non-attacking bishops.
    # Bishops on squares of the same color cannot attack each other.
    # The edge squares are evenly split between white and black.
    # To place the maximum number of bishops, we fill all edge squares of one color.
    max_bishops_placed = total_edge_squares // 2
    print(f"To place the maximum number of non-attacking bishops, we can place them on all edge squares of a single color.")
    print(f"The maximum number of bishops that can be placed is {max_bishops_placed}.")

    # Step 3: Calculate the number of empty edge squares.
    empty_edge_squares = total_edge_squares - max_bishops_placed
    print("\nThe number of edge squares that would lack bishops is calculated as:")
    print(f"{total_edge_squares} (total) - {max_bishops_placed} (bishops) = {empty_edge_squares}")

solve_bishops_puzzle()
<<<14>>>