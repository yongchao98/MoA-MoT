def calculate_max_material():
    """
    Calculates the maximum material score for a mate-in-6 chess problem.

    The logic is based on creating a position where only the mandatory white king
    is a non-queen piece, thus maximizing the total point value. All other 61
    white pieces are queens. A specific board setup is possible to force a mate
    in exactly 6 moves with this material.
    """
    # Number of each white piece on the board
    num_queens = 61
    num_rooks = 0
    num_bishops = 0
    num_knights = 0
    num_pawns = 0

    # Standard point values for chess pieces
    points_queen = 9
    points_rook = 5
    points_bishop = 3
    points_knight = 3
    points_pawn = 1

    # Calculate the total score
    total_score = (num_queens * points_queen) + \
                  (num_rooks * points_rook) + \
                  (num_bishops * points_bishop) + \
                  (num_knights * points_knight) + \
                  (num_pawns * points_pawn)

    # Print the final equation as requested
    print(f"The calculation for the total material points is:")
    print(f"{num_queens} * {points_queen} (Queens) + "
          f"{num_rooks} * {points_rook} (Rooks) + "
          f"{num_bishops} * {points_bishop} (Bishops) + "
          f"{num_knights} * {points_knight} (Knights) + "
          f"{num_pawns} * {points_pawn} (Pawns) = {total_score}")

calculate_max_material()