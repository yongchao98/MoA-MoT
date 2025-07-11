def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard after placing the maximum number of non-attacking bishops.
    """
    # An 8x8 board has 28 edge squares (8 on top + 8 on bottom + 6 on left + 6 on right).
    total_edge_squares = 28

    # To place the maximum number of non-attacking bishops, you must place them
    # on squares of the same color. The edge squares are evenly split between
    # white and black.
    # So, we can place a bishop on every white edge square, or every black one.
    max_bishops_on_edge = total_edge_squares / 2
    
    # The number of empty squares is the total minus the placed bishops.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    print(f"A standard 8x8 chessboard has {total_edge_squares} edge squares.")
    print(f"The maximum number of non-attacking bishops that can be placed on the edge is {int(max_bishops_on_edge)}.")
    print("This leaves the other squares empty.")
    print(f"Therefore, the number of edge squares that would lack bishops is:")
    print(f"{total_edge_squares} - {int(max_bishops_on_edge)} = {int(empty_edge_squares)}")

solve_bishops_puzzle()
<<<14>>>