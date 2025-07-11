def solve_bishops_on_edge():
    """
    Calculates how many edge squares on a chessboard would lack bishops
    if the maximum number of non-attacking bishops were placed on the edge.
    """

    # Step 1: Calculate the total number of edge squares on an 8x8 board.
    # An 8x8 board has 64 squares. The inner 6x6 grid has 36 squares.
    # So, the number of edge squares is 64 - 36 = 28.
    # Alternatively, it's 4 sides of 8 squares minus the 4 double-counted corners: 4*8 - 4 = 28.
    total_edge_squares = 28

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # A simple and optimal configuration is to use the first and last columns ('a' and 'h').
    # All 16 squares on these two columns are edge squares.
    # Bishops on the same column cannot attack each other.
    # We only need to check for attacks between column 'a' and column 'h'.
    # A bishop at (a1) attacks (h8). A bishop at (a8) attacks (h1).
    # These are the only two pairs of attacking bishops in this configuration.
    # To resolve these two conflicts, we must remove one bishop from each pair.
    # For example, we can remove the bishops at a1 and a8.
    # This leaves 16 - 2 = 14 non-attacking bishops. This is the maximum.
    max_bishops_on_edge = 14

    # Step 3: Calculate the number of empty edge squares.
    # This is the total number of edge squares minus the number of bishops placed.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    # Step 4: Print the final equation and the result.
    print("This problem asks for the number of empty edge squares when the maximum number of non-attacking bishops are placed.")
    print(f"Total number of edge squares on an 8x8 board: {total_edge_squares}")
    print(f"Maximum number of non-attacking bishops that can be placed on these edge squares: {max_bishops_on_edge}")
    print("\nThe number of edge squares that would lack bishops is the total minus the number of bishops placed:")
    print(f"{total_edge_squares} - {max_bishops_on_edge} = {empty_edge_squares}")

solve_bishops_on_edge()
<<<14>>>