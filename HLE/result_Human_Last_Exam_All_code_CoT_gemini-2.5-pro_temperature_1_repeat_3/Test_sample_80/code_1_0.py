def solve_bishop_puzzle():
    """
    Calculates the number of empty edge squares after placing the maximum
    number of non-attacking bishops on the edge of a chessboard.
    """
    board_size = 8

    # Step 1: Calculate the total number of edge squares.
    # This is the total number of squares minus the inner squares.
    total_squares = board_size * board_size
    inner_squares = (board_size - 2) * (board_size - 2)
    edge_squares = total_squares - inner_squares

    # Step 2: Determine the maximum number of non-attacking bishops
    # that can be placed on the edge squares.
    # A maximal placement can be achieved by placing bishops on one full edge (e.g., column 'a')
    # and the opposite edge excluding the corners (e.g., column 'h' except h1 and h8).
    # This gives N + (N-2) bishops.
    max_bishops = board_size + (board_size - 2)

    # Step 3: Calculate how many edge squares are left without bishops.
    empty_squares = edge_squares - max_bishops

    # Step 4: Print the explanation and the final equation.
    print("This puzzle asks for the number of empty edge squares on a chessboard after placing the maximum possible number of non-attacking bishops on them.")
    print(f"\n1. First, we find the total number of edge squares on an {board_size}x{board_size} board: {edge_squares}")
    print(f"2. Next, we find the maximum number of non-attacking bishops we can place on those squares: {max_bishops}")
    print("\n3. Finally, we find the number of squares that lack bishops by subtracting the number of bishops from the total number of edge squares.")
    print("\nThe final equation is:")
    print(f"{edge_squares} - {max_bishops} = {empty_squares}")

solve_bishop_puzzle()
<<<14>>>