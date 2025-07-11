def solve_chess_puzzle():
    """
    This function prints the move sequence for the chess puzzle.
    """
    # The best sequence of moves for black to force checkmate.
    # This pattern is known as Philidor's Legacy or a smothered mate.
    # 1... Black Knight checks the king, forcing it to move.
    # 2... The second Knight delivers a double check, forcing the king back.
    # 3... Black sacrifices the Queen, forcing the Rook to capture.
    # 4... The Rook's capture traps its own king, allowing the Knight to deliver the final checkmate.
    mate_sequence = "1... Nf2+ 2.Kg1 Nh3+ 3.Kh1 Qg1+ 4.Rxg1 Nf2#"
    print(mate_sequence)

solve_chess_puzzle()