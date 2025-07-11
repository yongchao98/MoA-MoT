def solve_bishops_puzzle():
    """
    Calculates how many edge squares would lack bishops if the maximum
    number of non-attacking bishops are placed on them.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # For an N x N board, this is 4 * (N - 1).
    total_edge_squares = 4 * (board_size - 1)

    # Step 2: Divide the squares by color.
    # On a standard chessboard, edge squares are split evenly.
    white_edge_squares = total_edge_squares // 2
    black_edge_squares = total_edge_squares // 2

    # Step 3: Analyze placement on white squares.
    # The white corners a8 and h1 are on the same diagonal.
    # This creates one conflict, so one square must remain empty.
    white_conflicts = 1
    max_white_bishops = white_edge_squares - white_conflicts

    # Step 4: Analyze placement on black squares.
    # The black corners a1 and h8 are on the same diagonal.
    # This creates one conflict, so one square must remain empty.
    black_conflicts = 1
    max_black_bishops = black_edge_squares - black_conflicts

    # Step 5: Calculate the total bishops and empty squares.
    total_bishops_placed = max_white_bishops + max_black_bishops
    empty_edge_squares = total_edge_squares - total_bishops_placed

    # Step 6: Print the explanation and the final calculation.
    print(f"A chessboard has {total_edge_squares} edge squares ({white_edge_squares} white and {black_edge_squares} black).")
    print("To place the maximum number of non-attacking bishops, we analyze each color separately.")
    print(f"For the {white_edge_squares} white edge squares, one bishop must be omitted because two squares ('a8' and 'h1') are on the same diagonal.")
    print(f"For the {black_edge_squares} black edge squares, one bishop must be omitted because two squares ('a1' and 'h8') are on the same diagonal.")
    print(f"\nMaximum bishops that can be placed: ({white_edge_squares} - {white_conflicts}) + ({black_edge_squares} - {black_conflicts}) = {total_bishops_placed}")
    print("\nThe number of edge squares lacking bishops is the total number of edge squares minus the total bishops placed.")
    print(f"Final Equation: {total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}")

solve_bishops_puzzle()
<<<2>>>