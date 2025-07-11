def solve_chess_puzzle():
    """
    This function prints the step-by-step solution to the chess puzzle.
    The solution is a forced checkmate sequence.
    """
    winning_sequence = [
        "1. Qxf7+ Kh7",
        "2. Bg6+ Kh8",
        "3. Nf7+ Kg8",
        "4. Nxh6++ Kh8",
        "5. Qg8+ Rxg8",
        "6. Nf7#"
    ]

    print("The winning sequence is:")
    for move in winning_sequence:
        print(move)

    # The number of moves for White to win is the number of moves White makes.
    white_moves_count = len(winning_sequence)
    print(f"\nWhite wins in {white_moves_count} moves.")

solve_chess_puzzle()