def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    total_squares = board_size * board_size
    inner_board_size = board_size - 2
    inner_squares = inner_board_size * inner_board_size
    total_edge_squares = total_squares - inner_squares

    # Step 2: Note the number of squares per color (not strictly needed for the calculation
    # but good for understanding). Total edge squares are split evenly.
    white_edge_squares = total_edge_squares // 2
    black_edge_squares = total_edge_squares // 2

    # Steps 3 & 4: The maximum number of non-attacking bishops on the
    # 14 white edge squares is 7. By symmetry, it's also 7 for the 14 black edge squares.
    max_bishops_on_white_edges = 7
    max_bishops_on_black_edges = 7

    # Step 5: Calculate the total number of bishops that can be placed.
    total_bishops_placed = max_bishops_on_white_edges + max_bishops_on_black_edges

    # Step 6: Calculate the number of edge squares lacking bishops.
    empty_edge_squares = total_edge_squares - total_bishops_placed

    print("Problem: How many edge squares would lack bishops if the maximum number of non-attacking bishops are placed on them?")
    print(f"\n1. Total edge squares on a {board_size}x{board_size} board: {total_edge_squares}")
    print(f"2. Maximum non-attacking bishops that can be placed on white edge squares: {max_bishops_on_white_edges}")
    print(f"3. Maximum non-attacking bishops that can be placed on black edge squares: {max_bishops_on_black_edges}")
    print(f"4. Total bishops placed = {max_bishops_on_white_edges} + {max_bishops_on_black_edges} = {total_bishops_placed}")
    print("\nFinal Calculation:")
    print(f"The number of edge squares that would lack bishops is the total number of edge squares minus the number of bishops placed.")
    print(f"Result: {total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}")

solve_chess_puzzle()
<<<14>>>