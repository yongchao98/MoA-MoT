def solve_chess_puzzle():
    """
    This function analyzes the given chess position and determines the best move for White.
    """
    position_fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    best_move = "Qxc5+"
    
    explanation = """
    The best move for White is to capture the black queen.
    
    1.  The move `1. d7` seems tempting to promote the pawn, but it fails to `1... Qc1+`, leading to a forced queen trade and a losing pawn endgame for White.
    2.  The correct move is `1. Qxc5+`.
    3.  This move is a check, and since Black's king on b6 cannot recapture, Black is forced to move the king to the side of the board (either `Ka6` or `Kb7`).
    4.  This results in a Queen vs Pawns endgame that is completely winning for White. White's queen dominates the board, supports the advance of the d-pawn, and targets the weak black king.
    """
    
    print("Analysis of the position:")
    print(explanation)
    print("The best move for White is:")
    print(best_move)

solve_chess_puzzle()