def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_side_length = 8

    # Step 1: Calculate the total number of edge squares.
    # An 8x8 board has 4 sides. If we simply multiply 8 * 4, we count
    # the 4 corner squares twice. So, we subtract 4.
    total_edge_squares = (board_side_length * 4) - 4

    # Step 2: Determine the maximum number of non-attacking bishops.
    # To place the maximum number of non-attacking bishops, they must all
    # be on squares of the same color. The 28 edge squares of a chessboard
    # are composed of 14 light and 14 dark squares.
    # Therefore, we can place a maximum of 14 bishops on the edge.
    max_bishops_on_edge = total_edge_squares // 2

    # Step 3: Calculate the number of empty edge squares.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    # Step 4: Print the final equation with all numbers.
    print("This problem can be solved with the following calculation:")
    print(f"Total edge squares on an 8x8 board: ({board_side_length} * 4) - 4 = {total_edge_squares}")
    print(f"Maximum non-attacking bishops that can be placed on the edge: {total_edge_squares} / 2 = {max_bishops_on_edge}")
    print("The number of squares lacking bishops is the total edge squares minus the placed bishops.")
    print(f"Final equation: {total_edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")

solve_bishop_puzzle()
<<<14>>>