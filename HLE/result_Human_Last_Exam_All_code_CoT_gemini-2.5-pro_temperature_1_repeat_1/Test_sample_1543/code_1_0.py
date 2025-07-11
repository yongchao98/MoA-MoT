def solve_chess_puzzle():
    """
    This function determines the minimal number of moves for White to win.
    
    The position is given by the FEN: 9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1
    
    1. White plays Qe4.
       - This move controls the black king's only escape square, i8.
       - It threatens mate in one with Qj4#.
    
    2. Black must respond. We analyze the main defensive options:
       - If Black plays 1... cf5 (blocking the queen's control of i8), White responds with 2. Ai4#, which is checkmate.
         The Archbishop (knight move) checks the king on j8. The king cannot escape to i8 or j7 (covered by the Queen on e4). The check is unblockable and the Archbishop cannot be captured.
       - If Black plays any other move that does not stop the threat (e.g. 1... ph5), White plays 2. Qj4#, which is checkmate.
    
    Therefore, White can force a checkmate in 2 moves.
    """
    
    # The minimal number of moves for White to win.
    moves_to_win = 2
    
    print(moves_to_win)

solve_chess_puzzle()