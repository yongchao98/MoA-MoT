def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing
    the maximum possible number of non-attacking bishops on the edge.
    """

    # Step 1: Calculate the total number of edge squares on an 8x8 board.
    # An 8x8 board has 8 squares on the top and bottom rows, and 6 unique
    # squares on the left and right columns (excluding corners).
    total_edge_squares = 8 + 8 + 6 + 6

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # Bishops attack along diagonals of the same color.
    # The 28 edge squares are split evenly: 14 light and 14 dark.
    # To place the maximum number of non-attacking bishops, we can fill all
    # edge squares of a single color.
    max_bishops_placed = total_edge_squares // 2

    # Step 3: Calculate how many edge squares are left without bishops.
    empty_edge_squares = total_edge_squares - max_bishops_placed

    # Step 4: Print the final equation as requested.
    print(f"To find the number of empty edge squares, we subtract the maximum number of bishops placed from the total number of edge squares.")
    print(f"The final calculation is:")
    print(f"{total_edge_squares} - {max_bishops_placed} = {empty_edge_squares}")


solve_chess_puzzle()
<<<14>>>