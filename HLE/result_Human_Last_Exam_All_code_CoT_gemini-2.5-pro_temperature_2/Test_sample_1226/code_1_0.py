def calculate_max_material():
    """
    Calculates the maximum material points for a constructed chess position
    where White can mate in exactly 6 moves, but not fewer.
    The calculation is based on a path-clearing mechanism where the board
    is filled with pieces to prevent any faster mates.
    """
    # Standard piece point values
    queen_points = 9
    rook_points = 5
    bishop_points = 3
    knight_points = 3
    pawn_points = 1

    # Total squares on a chessboard
    total_squares = 64
    # Squares occupied by the two kings
    king_squares = 2
    # Minimum required empty squares for the solution to work:
    # 1 square for the black king to shuttle to.
    # 5 squares for the blocker pieces to move to.
    required_empty_squares = 6

    # Number of squares available for white pieces
    white_piece_squares = total_squares - king_squares - required_empty_squares

    # Define the pieces essential for the mate-in-6 mechanism.
    # These pieces cannot be Queens (except the mating one) due to their specific roles.
    num_pawns = 2  # Part of the cage, must not attack king's starting rank
    num_knights = 1  # Part of the cage, covers key escape squares

    # The remaining pieces can all be Queens to maximize points.
    # This includes the mating queen, the path-blocking queens, and the filler queens.
    num_rooks = 0
    num_bishops = 0
    num_queens = white_piece_squares - num_pawns - num_knights
    
    # Calculate the total score
    total_points = (num_queens * queen_points) + \
                   (num_rooks * rook_points) + \
                   (num_bishops * bishop_points) + \
                   (num_knights * knight_points) + \
                   (num_pawns * pawn_points)

    print("Calculation of maximum material points for a mate in 6:")
    print("-" * 60)
    print(f"{num_queens} queens * {queen_points} points = {num_queens * queen_points}")
    print(f"{num_rooks} rooks * {rook_points} points = {num_rooks * rook_points}")
    print(f"{num_bishops} bishops * {bishop_points} points = {num_bishops * bishop_points}")
    print(f"{num_knights} knight * {knight_points} points = {num_knights * knight_points}")
    print(f"{num_pawns} pawns * {pawn_points} points = {num_pawns * pawn_points}")
    print("-" * 60)
    print(f"Total points = {total_points}")


calculate_max_material()