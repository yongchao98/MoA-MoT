def solve_chess_puzzle():
    """
    Analyzes the chess position and provides the best move for White.
    FEN: 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    """
    print("Analyzing the chess position to find the best move for White.")
    print("=" * 30)

    print("Initial assessment:")
    print("White's strength is the passed d6-pawn. Black's queen and king are defending.")
    print("A quiet move fails to ...Qxd6, so White must play a check.")
    print("\n")

    print("Evaluating the trap move 1. Qxc5+:")
    print("1. Qxc5+ bxc5 leads to a King and Pawn endgame.")
    print("Although it seems simple, Black's king is close enough to stop the d-pawn, and the new c-pawn gives Black a winning advantage.")
    print("Therefore, 1. Qxc5+ is a losing move.")
    print("\n")
    
    print("The correct strategy is to force the Black king away from the defense.")
    print("The winning move is 1. Qb2+.")
    print("\n")

    print("Here is the main winning line:")
    print("Move 1: White plays Qb2+")
    print("White's move: 1. Qb2+")
    print("This check forces the Black king to the a-file, away from the key c7 square.")
    print("Black's reply: 1... Ka6 (or 1... Ka7)")
    print("-" * 20)

    print("Move 2: White plays d7")
    print("White's move: 2. d7")
    print("This advances the pawn, forcing Black's queen to pin it.")
    print("Black's reply: 2... Qe7")
    print("-" * 20)
    
    print("Move 3: White plays Qb5+")
    print("White's move: 3. Qb5+")
    print("This is the decisive blow. It forces the Black king to a7 (as ...Ka5 would be mate).")
    print("Black's reply: 3... Ka7")
    print("-" * 20)

    print("\nFinal Position Analysis:")
    print("After this sequence, Black is paralyzed. The king is stuck on the a-file, and the queen is tied to the defense of the d7-pawn.")
    print("White will simply march the king over (Kg2-f3-e4...) to break the pin and ensure the pawn's promotion.")
    print("Black has no counterplay and cannot stop this plan.")

    best_move = "Qb2+"
    print(f"\nConclusion: White's best move is {best_move}.")

solve_chess_puzzle()
<<<Qb2+>>>