def solve_bishops_puzzle():
    """
    Calculates the number of empty edge squares on a chessboard when the
    maximum number of non-attacking bishops are placed on the edge.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # An N x N board has 4 sides of N squares, with 4 corners counted twice.
    # Formula: (N * 4) - 4
    total_edge_squares = (board_size * 4) - 4

    # Step 2: Determine the maximum number of non-attacking bishops on the edge.
    # The edge squares are divided equally between two colors (light and dark).
    # On an 8x8 board, there are 14 light and 14 dark edge squares.
    # Bishops can't attack bishops on squares of a different color.
    # So, the maximum number of non-attacking bishops we can place is the
    # number of edge squares of a single color.
    max_bishops_on_edge = total_edge_squares / 2

    # Step 3: Calculate how many edge squares are empty.
    # This is the total number of edge squares minus the ones occupied by bishops.
    empty_edge_squares = total_edge_squares - max_bishops_on_edge

    print("To solve this, we first find the total number of edge squares on an 8x8 board.")
    print(f"An 8x8 board has {total_edge_squares} edge squares.")
    print("\nTo place the maximum number of non-attacking bishops, we can occupy all edge squares of a single color.")
    print(f"There are {int(max_bishops_on_edge)} edge squares of one color.")
    print(f"Therefore, a maximum of {int(max_bishops_on_edge)} bishops can be placed on the edge without attacking each other.")
    print("\nThe number of edge squares that lack bishops is the total minus the number of bishops placed.")
    print("\nFinal calculation:")
    print(f"{total_edge_squares} - {int(max_bishops_on_edge)} = {int(empty_edge_squares)}")

solve_bishops_puzzle()

# The final answer is the number of empty squares.
print("\n<<<14>>>")