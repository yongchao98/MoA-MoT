def solve_chess_puzzle():
    """
    This function prints the winning chess move and variation.
    The position is a tense endgame where White must be precise to convert the advantage of the a-pawn.
    The best move is to immediately promote the pawn, sacrificing it for the defending knight.
    This simplifies the position into a winning Knight vs Pawns endgame.
    """
    move_1_white = "a8=Q"
    move_1_black = "Nxa8"
    move_2_white = "Nc5"

    # The equation represents the move sequence.
    # 1. White promotes the pawn, forcing Black to capture.
    # 2. White then repositions the knight to a dominant square (c5),
    #    attacking a pawn and preparing to neutralize Black's counterplay.
    #    This leads to a winning position for White.
    print(f"1. {move_1_white} {move_1_black} 2. {move_2_white}")

solve_chess_puzzle()