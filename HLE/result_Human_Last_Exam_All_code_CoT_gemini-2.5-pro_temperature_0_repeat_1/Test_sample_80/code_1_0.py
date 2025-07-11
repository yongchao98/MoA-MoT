def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # Each of the 4 sides has (board_size - 1) unique edge squares.
    total_edge_squares = 4 * (board_size - 1)
    print(f"An {board_size}x{board_size} chessboard has {total_edge_squares} edge squares.")

    # Step 2: The problem can be split by color.
    # The edge squares are evenly divided between white and black.
    edge_squares_per_color = total_edge_squares // 2
    print(f"These are divided into {edge_squares_per_color} white squares and {edge_squares_per_color} black squares.")

    # Step 3: Determine the maximum number of non-attacking bishops on one color's edge squares.
    # Through analysis, it's found that you can place a maximum of 7 non-attacking
    # bishops on the 14 white edge squares (e.g., on a1, c1, e1, g1, b8, d8, f8).
    max_bishops_one_color = 7
    print(f"On the {edge_squares_per_color} squares of a single color, a maximum of {max_bishops_one_color} non-attacking bishops can be placed.")

    # Step 4: Calculate the total number of bishops that can be placed.
    # This is the sum for both white and black squares.
    total_bishops_placed = max_bishops_one_color * 2
    print(f"By symmetry, the total number of bishops that can be placed on all edge squares is {max_bishops_one_color} + {max_bishops_one_color} = {total_bishops_placed}.")

    # Step 5: Calculate the number of empty edge squares.
    empty_edge_squares = total_edge_squares - total_bishops_placed
    print("\nThe number of edge squares lacking bishops is the total number of edge squares minus the number of bishops placed.")
    print(f"Final Equation: {total_edge_squares} - {total_bishops_placed} = {empty_edge_squares}")

solve_bishops_puzzle()
<<<14>>>