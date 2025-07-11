def solve_chess_puzzle():
    """
    This function prints the mate-in-2 sequence for the given chess puzzle.
    
    The logic is as follows:
    1. Black's first move is Knight to e3 (Ne3). This is not a check but creates an
       unstoppable threat.
    2. If White does not capture the knight, Black's next move is Knight to f1,
       which is checkmate (Nf1#). This is the main line.
    3. If White's only defense is to capture the knight (Bxe3), Black's response is
       Rook to f1, which is also checkmate (Rf1#).
    4. Since a mate is guaranteed on the second move regardless of White's play,
       the puzzle is solved. The printed sequence shows the primary threat.
    """
    first_move = "Ne3"
    second_move = "Nf1#"
    
    print(f"{first_move}, {second_move}")

solve_chess_puzzle()
