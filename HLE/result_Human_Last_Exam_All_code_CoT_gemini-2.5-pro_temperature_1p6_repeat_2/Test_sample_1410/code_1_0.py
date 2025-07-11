def solve_chess_puzzle():
    """
    This function prints the move-by-move solution to the chess puzzle.
    The FEN for the position is: 4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20
    It's a mate in 5 for White, assuming optimal play from Black.
    """
    
    # White plays optimally to force a checkmate.
    # Black plays optimally to delay the checkmate as long as possible.
    # The solution is a mate in 5 moves.

    moves = [
        ("1.", "Qxf7+", "Kh8"),
        ("2.", "Nxg6+", "fxg6"),
        ("3.", "Bxg6+", "Kg8"),
        ("4.", "Qe8+", "Rxe8"),
        ("5.", "Rf8#", "")
    ]
    
    print("The solution is a forced mate in 5 moves for White.")
    print("The sequence, assuming optimal play from Black, is:")
    print("-" * 30)
    print("Move |  White  |  Black")
    print("-" * 30)
    for move in moves:
        move_num, white_move, black_move = move
        print(f"{move_num:<4} | {white_move:<7} | {black_move}")
    print("-" * 30)


solve_chess_puzzle()