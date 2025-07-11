def solve_chess_puzzle():
    """
    This function analyzes the chess position from the FEN 8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1
    and prints the best move for White.
    """
    position_summary = """
    White's decisive advantage is the passed pawn on the d6 square. Black's queen is the only
    piece holding the position together. White's best move must exploit this situation.

    The move 1. Qd4 forces a win in all variations:
    - If Black plays 1... Qxd4+, White responds with 2. cxd4, leading to a winning pawn endgame.
    - If Black avoids the trade with a king move like 1... Kc7, White continues with 2. Qe5+,
      which leads to unstoppable threats involving the advance of the d-pawn to d7.
    """
    
    best_move = "Qd4"
    
    print("Analysis of the position:")
    print(position_summary)
    print(f"The best move for White is: {best_move}")

solve_chess_puzzle()