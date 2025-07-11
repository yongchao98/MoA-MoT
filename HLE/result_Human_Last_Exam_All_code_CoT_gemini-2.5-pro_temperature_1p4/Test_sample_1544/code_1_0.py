def solve_chess_puzzle():
    """
    This function analyzes the given chess position and prints the step-by-step reasoning for White's best move.
    """
    print("--- Chess Puzzle Analysis ---")
    print("Position FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1")
    print("")

    print("Step 1: Understanding the Position")
    print("White has a major advantage with a passed pawn on d6, but it is attacked by Black's active queen on c5.")
    print("White must find a forcing move to convert this advantage into a win.")
    print("")

    print("Step 2: Evaluating a plausible but incorrect move")
    print("The move 1. Qxc5+ seems logical to trade off Black's dangerous queen. However, after 1... bxc5, the resulting pawn endgame is losing for White.")
    print("Black's king is too close and can easily stop the d-pawn. For example: 2. d7 Kc7 3. d8=Q+ Kxd8, and Black wins.")
    print("")

    print("Step 3: Identifying the winning move")
    print("The correct move is 1. Qf3. This attacks the Black queen, forcing a favorable sequence for White.")
    print("Black's best reply is to trade queens, leading to a King and Pawn endgame that White can win.")
    print("")

    print("Step 4: The Winning Variation")
    print("The game demonstrates White's winning plan after 1. Qf3:")
    print("1. Qf3 Qxf3+")
    print("2. Kxf3 Kc5   (Black's best attempt, activating the king to attack the pawn)")
    print("3. Ke4   Kb4")
    print("4. Kd4   Kxa4")
    print("5. Kc5!  (This is the decisive move. It places Black's king in zugzwang, forcing it away from the key defensive squares.)")
    print("\nAfter this, White's d-pawn will advance to promotion, securing the win.")
    print("")
    
    winning_move = "Qf3"
    print(f"Conclusion: The best move for White is {winning_move}.")

solve_chess_puzzle()
<<<Qf3>>>