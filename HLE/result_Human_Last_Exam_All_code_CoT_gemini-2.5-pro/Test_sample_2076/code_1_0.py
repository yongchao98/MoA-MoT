def solve_chess_puzzle():
    """
    This function determines the minimal pieces for a Diagonal Corridor Mate
    with the given constraints and prints the result in the specified format.
    """

    # The problem requires constructing a Diagonal Corridor Mate.
    # The analysis leads to a minimal 4-piece solution:
    # 1. White Bishop (checking piece)
    # 2. White Pawn (protecting the Bishop)
    # 3. Black Pawn (blocking an escape square)
    # 4. Black Knight (blocking the other escape square, cannot be a pawn on the 8th rank)
    # The Knight and Bishop have the same value (3), so either is valid for the minimal
    # value solution. We choose the Knight.

    # Define the required pieces for each color.
    white_pieces = ["Pawn", "Bishop"]
    black_pieces = ["Pawn", "Knight"]

    # The problem asks for sorting by piece value.
    # Piece values: Pawn=1, Knight=3, Bishop=3, Rook=5, Queen=9.
    # The lists are already sorted according to these values.

    # Format the lists by adding the color prefix.
    formatted_white_pieces = [f"White {piece}" for piece in white_pieces]
    formatted_black_pieces = [f"Black {piece}" for piece in black_pieces]

    # Combine the lists, with White's pieces first.
    final_piece_list = formatted_white_pieces + formatted_black_pieces

    # Join the list into a single, comma-separated string.
    result_string = ", ".join(final_piece_list)

    print(result_string)

solve_chess_puzzle()