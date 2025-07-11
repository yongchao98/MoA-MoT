def solve_capablanca_chess_mate():
    """
    This function analyzes the given Capablanca chess position and determines the minimum moves to win.
    The logic is explained in the text accompanying the code.
    """

    fen = "9k/5c1pb1/10/10/10/3Q6/PP5A2/K9 w - - 0 1"
    
    # Analysis of the position leads to a mate in 2.
    # The line is:
    # 1. Qj3+  (White Queen from d3 to j3).
    #    This is a double attack:
    #    a) It checks the Black King on j8.
    #    b) It attacks the Black Chancellor on f7.
    #    Black is forced to move the King. The best response is Kh8.
    #    (Ki8 is also possible but leads to a similar mate).
    #
    # 2. Qxf7# (White Queen captures Chancellor on f7, delivering checkmate).
    #    a) The Queen on f7 now checks the King on h8 diagonally.
    #    b) The King's escape squares (g8, g7) are both controlled by the Queen on f7.
    #    c) The Queen cannot be captured, and the check cannot be blocked.
    #
    # The game ends in 2 moves for White.

    moves_to_win = 2
    
    print("The chess position is described by the FEN: {}".format(fen))
    print("White must make a move, then Black responds, then White delivers the final blow.")
    print("The winning sequence is as follows:")
    print("1. White moves Queen from d3 to j3 (check and attack).")
    print("1. ... Black is forced to move King from j8 to h8.")
    print("2. White moves Queen to capture the Chancellor on f7, which is checkmate.")
    print("\nThis is a mate in 2 moves.")
    print("Minimal amount of moves by White to win: {}".format(moves_to_win))


solve_capablanca_chess_mate()
