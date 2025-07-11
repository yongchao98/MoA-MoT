def solve_bishops_puzzle():
    """
    Calculates how many edge squares lack bishops when the maximum number of
    non-attacking bishops are placed on the edge of an 8x8 chessboard.
    """
    # Step 1: Calculate the total number of edge squares.
    # An 8x8 board has 8 squares on the top row, 8 on the bottom row,
    # and 6 unique squares on each of the side columns (the corners are already counted).
    squares_top_row = 8
    squares_bottom_row = 8
    squares_left_side_inner = 6
    squares_right_side_inner = 6
    total_edge_squares = squares_top_row + squares_bottom_row + squares_left_side_inner + squares_right_side_inner

    # Step 2: Determine the maximum number of non-attacking bishops.
    # Bishops can only attack pieces on the same color square. The 28 edge squares
    # are split into 14 light squares and 14 dark squares. We can maximize
    # the bishops on each color set independently.

    # Through analysis of the attack patterns, the maximum number of non-attacking
    # bishops that can be placed on the 14 light-colored edge squares is 7.
    max_bishops_light_squares = 7

    # Similarly, the maximum number of non-attacking bishops that can be placed
    # on the 14 dark-colored edge squares is 7.
    max_bishops_dark_squares = 7

    # The total maximum number of bishops is the sum from both color sets.
    max_placed_bishops = max_bishops_light_squares + max_bishops_dark_squares

    # Step 3: Calculate the number of edge squares without bishops.
    empty_edge_squares = total_edge_squares - max_placed_bishops

    # Step 4: Print the explanation and the final equation.
    print(f"The total number of edge squares on an 8x8 board is {total_edge_squares}.")
    print(f"The maximum number of non-attacking bishops that can be placed on these squares is {max_placed_bishops}.")
    print("The number of squares lacking bishops is the total number of edge squares minus the number of placed bishops.")
    print("\nFinal Equation:")
    print(f"{total_edge_squares} - {max_placed_bishops} = {empty_edge_squares}")

solve_bishops_puzzle()
<<<14>>>