def solve_chess_puzzle():
    """
    Analyzes the given chess FEN and determines the best move for White.
    """
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"

    print("Analyzing the chess position with FEN:", fen)
    print("\n--- Step 1: Position Evaluation ---")
    print("White is in a difficult position. Black is up one pawn, and White's king is potentially vulnerable.")
    print("If White makes a quiet move, Black can play ...Qxf2+, forcing a queen trade.")
    print("This leads to a king and pawn endgame that is lost for White, as Black's king on b6 is far more active and can attack White's pawns.")

    print("\n--- Step 2: Analyzing a Losing Candidate Move ---")
    print("Let's consider the immediate queen trade: 1. Qxc5+")
    print("Black must recapture with the king: 1... Kxc5")
    print("This results in a pawn endgame where Black's king on c5 is perfectly placed to win White's pawns. White's king on g1 is too far away to defend. Therefore, 1. Qxc5+ leads to a loss.")

    print("\n--- Step 3: Finding the Best Strategy (Forcing a Draw) ---")
    print("Since passive moves and trading queens both lose, White's only hope is to find a forcing sequence of checks that never allows Black to escape.")
    print("This is known as a perpetual check, and it results in a draw.")
    print("The most precise move to start this sequence is 1. Qf8+.")

    print("\n--- Step 4: The Winning (Drawing) Sequence ---")
    print("The forced sequence is as follows:")
    print("Move 1: White plays Qf8+. Black's king is forced to move to a6 (1... Ka6).")
    print("Move 2: White plays Qa8+. Black's king is forced to move back to b6 (2... Kb6).")
    print("Move 3: White can now play Qf8+ again, repeating the position. Black cannot avoid the checks, so the game is a draw by repetition.")
    
    print("\n--- Conclusion ---")
    print("Securing a draw is the best possible outcome for White.")
    final_move = "Qf8+"
    print(f"The best move for White is: {final_move}")

solve_chess_puzzle()