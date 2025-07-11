def solve_bishops_on_edge():
    """
    Calculates the number of empty edge squares on a chessboard after placing
    the maximum number of non-attacking bishops on the edge.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An n x n board has 4 sides of n squares, but the 4 corners are counted twice.
    # So, the formula is (n * 4) - 4.
    total_edge_squares = (board_size * 4) - 4

    # Step 2: Determine the maximum number of non-attacking bishops.
    # Bishops on white squares do not attack bishops on black squares.
    # The edge squares of a chessboard are split evenly by color.
    # To place the maximum number of bishops, we fill all edge squares of one color.
    max_bishops_on_edge = total_edge_squares // 2

    # Step 3: Calculate how many edge squares are empty.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    # Step 4: Print the reasoning and the final equation.
    print(f"First, we find the total number of edge squares on an {board_size}x{board_size} board.")
    print(f"Total edge squares = ({board_size} * 4) - 4 = {total_edge_squares}")
    print(f"To place the maximum number of non-attacking bishops, we can place them on all edge squares of a single color.")
    print(f"Maximum bishops that can be placed on the edge = {total_edge_squares} / 2 = {max_bishops_on_edge}")
    print(f"The number of edge squares lacking bishops is the total minus the number placed.")
    print(f"Final calculation: {total_edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")

solve_bishops_on_edge()
<<<14>>>