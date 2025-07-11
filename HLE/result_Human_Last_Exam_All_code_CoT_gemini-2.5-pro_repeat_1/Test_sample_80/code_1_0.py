def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of an 8x8 chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An N x N board has 4 sides of N squares, but the 4 corners are counted twice.
    # So, the formula is 4 * N - 4.
    total_edge_squares = 4 * board_size - 4

    # Step 2: The edge squares are split evenly between two colors.
    white_edge_squares = total_edge_squares // 2
    black_edge_squares = total_edge_squares // 2

    # Step 3: For each color, find the maximum number of non-attacking bishops.
    # On the white edge squares, only the corners a1 and h8 attack each other.
    # So we must leave one of them empty.
    # Number of conflicting pairs on white squares is 1.
    max_white_bishops = white_edge_squares - 1

    # On the black edge squares, only the corners a8 and h1 attack each other.
    # So we must leave one of them empty.
    # Number of conflicting pairs on black squares is 1.
    max_black_bishops = black_edge_squares - 1

    # Step 4: The total maximum number of bishops is the sum for both colors.
    total_max_bishops = max_white_bishops + max_black_bishops

    # Step 5: The number of empty squares is the total minus the number of placed bishops.
    empty_edge_squares = total_edge_squares - total_max_bishops

    print("This puzzle is about placing the maximum number of non-attacking bishops on the edge of a chessboard.")
    print(f"A standard {board_size}x{board_size} board has {total_edge_squares} edge squares.")
    print(f"The maximum number of non-attacking bishops that can be placed on these squares is {total_max_bishops}.")
    print("The number of edge squares that would lack bishops is the total number of edge squares minus the number of placed bishops.")
    print("\nThe final calculation is:")
    # The final equation as requested by the user
    print(f"{total_edge_squares} - {total_max_bishops} = {empty_edge_squares}")

solve_bishops_puzzle()
<<<2>>>