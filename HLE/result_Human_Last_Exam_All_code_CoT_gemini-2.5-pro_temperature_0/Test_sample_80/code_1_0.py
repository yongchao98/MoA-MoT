def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edges of a chessboard.
    """
    # Step 1: Calculate the total number of edge squares on an 8x8 board.
    # An 8x8 board has 4 sides, each with 8 squares. The 4 corners are counted twice.
    # Total edge squares = (8 * 4) - 4
    total_edge_squares = (8 * 4) - 4

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # As explained in the plan, we can place 8 bishops on the 'a' file and
    # 6 bishops on the 'h' file (from h2 to h7).
    bishops_on_a_file = 8
    bishops_on_h_file = 6
    max_bishops_on_edge = bishops_on_a_file + bishops_on_h_file

    # Step 3: Calculate the number of edge squares without bishops.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    # Print the final equation and the result.
    print("Total edge squares on a chessboard: 28")
    print("Maximum non-attacking bishops that can be placed on the edge: 14")
    print("The number of edge squares lacking bishops is calculated as:")
    print(f"{total_edge_squares} (total edge squares) - {max_bishops_on_edge} (placed bishops) = {empty_edge_squares}")

solve_bishops_puzzle()