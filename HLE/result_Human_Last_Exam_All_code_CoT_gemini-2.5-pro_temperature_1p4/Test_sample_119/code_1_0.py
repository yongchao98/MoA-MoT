def solve_chess_puzzle():
    """
    This function identifies and prints the shortest mating sequence for Black.
    """
    # The sequence is a mate-in-2, starting with Black's move.
    # 1. Black checks the White King.
    # 2. White has only one forced move in response.
    # 3. Black delivers the final checkmate.
    mating_sequence = ["h5+", "Kxh5", "f6#"]

    # Print the moves in standard chess notation, separated by spaces.
    print(' '.join(mating_sequence))

solve_chess_puzzle()