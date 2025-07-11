def solve_chess_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing
    the maximum possible number of non-attacking bishops on the edge.
    """
    # Step 1: Calculate the total number of edge squares on an 8x8 board.
    # There are 8 squares on each of the 4 sides, but the 4 corners are counted twice.
    side_length = 8
    corners = 4
    total_edge_squares = (side_length * 4) - corners

    # Step 2: Determine the maximum number of bishops that can be placed on the edge.
    # To maximize the number of non-attacking bishops, we can place them on all
    # edge squares of a single color. The 28 edge squares are split evenly
    # into 14 white and 14 black squares.
    max_bishops_placed = total_edge_squares // 2

    # Step 3: Calculate the number of edge squares that will be empty.
    # This is the total number of edge squares minus the number of squares occupied by bishops.
    empty_edge_squares = total_edge_squares - max_bishops_placed

    # Step 4: Print the explanation and the final equation.
    print("A standard 8x8 chessboard has 28 edge squares.")
    print("To place the maximum number of non-attacking bishops, you must occupy all edge squares of a single color.")
    print(f"There are {max_bishops_placed} white edge squares and {max_bishops_placed} black edge squares.")
    print("Therefore, the maximum number of bishops you can place on the edge is 14.")
    print("\nThe number of empty edge squares is the total number of edge squares minus the number of bishops placed.")
    print(f"Final calculation: {total_edge_squares} - {max_bishops_placed} = {empty_edge_squares}")

solve_chess_puzzle()
<<<14>>>