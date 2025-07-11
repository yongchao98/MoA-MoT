def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on them.
    """

    board_size = 8
    print(f"Answering for a standard {board_size}x{board_size} chessboard.")

    # Step 1: Calculate the total number of edge squares.
    # Total squares = 8*8 = 64. Inner squares = 6*6 = 36. Edge = 64 - 36 = 28.
    total_edge_squares = (board_size * 4) - 4
    print(f"The total number of edge squares is {total_edge_squares}.")

    # Step 2: Divide the problem by square color.
    # The edge squares are evenly split between light and dark colors.
    edge_squares_per_color = total_edge_squares // 2
    print(f"The edge is made of {edge_squares_per_color} light and {edge_squares_per_color} dark squares.")

    # Step 3: Determine the maximum number of non-attacking bishops per color set.
    # On the 14 light-colored edge squares, a maximal placement is 7 bishops
    # (e.g., b1, d1, f1, h1 on the top row, and c8, e8, g8 on the bottom row).
    # This configuration is non-attacking and covers all 14 light edge squares.
    max_bishops_on_one_color = 7
    print(f"The maximum number of non-attacking bishops that can be placed on the edge squares of a single color is {max_bishops_on_one_color}.")

    # Step 4: Calculate the total number of bishops that can be placed.
    # By symmetry, we can also place 7 bishops on the dark edge squares.
    total_bishops_placed = max_bishops_on_one_color + max_bishops_on_one_color
    print(f"So, the total maximum bishops placed is {max_bishops_on_one_color} + {max_bishops_on_one_color} = {total_bishops_placed}.")

    # Step 5: Find the number of edge squares without bishops.
    empty_edge_squares = total_edge_squares - total_bishops_placed

    print("\nThe final answer is the total edge squares minus the total bishops placed.")
    print(f"Result: The number of edge squares lacking bishops is {total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}.")

solve_bishops_puzzle()
<<<14>>>