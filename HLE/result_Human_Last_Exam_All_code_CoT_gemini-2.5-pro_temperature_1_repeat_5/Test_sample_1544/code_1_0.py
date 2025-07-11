def solve_chess_puzzle():
    """
    Analyzes the chess position from the FEN 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    and prints the reasoning to find the best move for White.
    """
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"

    print("Analyzing the chess position to find the best move for White.")
    print(f"FEN: {fen}\n")

    print("--- Step 1: Initial Position Analysis ---")
    print("White's key asset is the far-advanced passed pawn on d6.")
    print("However, this pawn is currently blocked by Black's pawn on d7.")
    print("Black's queen on c5 is extremely active and controls important squares, posing a threat to White.")
    print("White is down one pawn, so the position is dangerous and requires precise play.\n")

    print("--- Step 2: Evaluating Candidate Move 'Qxc5+' ---")
    print("One option is to trade queens with 1. Qxc5+.")
    print("After 1. Qxc5+, Black's king is in check. Black must recapture with the pawn on c6: 1...cxc5.")
    print("This trade leads to a king and pawn endgame.")
    print("The resulting pawn structure is:")
    print("  - White: Pawns on a4, d6")
    print("  - Black: Pawns on a5, c5, d7")
    print("Crucially, Black's king on b6 is very active, while White's king on g1 is far from the action.")
    print("Black's active king can attack and win White's a4 pawn, and then create a winning passed pawn on the queenside.")
    print("Conclusion: The move 1. Qxc5+ leads to a losing endgame for White.\n")

    print("--- Step 3: Evaluating Candidate Move 'Qxf7' ---")
    print("A second option is 1. Qxf7.")
    print("This move captures the d7-pawn, which unblocks White's own d6-pawn, creating an immediate promotion threat.")
    print("Black's best defense is to use the active queen to force a draw by giving endless checks (perpetual check).")
    print("A likely sequence of moves would be: 1...Qc1+ 2. Kg2 Qc2+.")
    print("Now, if White plays 3. Kh3, Black must continue checking with 3...Qc1+ to avoid losing.")
    print("Any other move by Black (like 3...Qxc6) would allow White's d-pawn to advance with 4. d7, leading to a win for White.")
    print("Thus, Black will force a draw by repetition.")
    print("Conclusion: The move 1. Qxf7 leads to a forced draw.\n")

    print("--- Step 4: Final Conclusion ---")
    print("We have analyzed two main moves:")
    print("  - 1. Qxc5+ leads to a certain loss.")
    print("  - 1. Qxf7 leads to a forced draw.")
    print("In chess, a draw is a much better result than a loss.")
    print("Therefore, playing for the draw is the optimal strategy for White.\n")

    best_move_san = "Qxf7"
    print(f"The best move for White is: {best_move_san}")

solve_chess_puzzle()