def solve_chess_mate():
    """
    This function prints the shortest mating sequence for Black from the described position.
    The sequence is a mate in 2.
    """
    # 1. Black sacrifices the h-pawn to force the White king onto the h-file.
    black_move_1 = "h6+"

    # 2. White has no other choice but to capture the pawn.
    white_move_1 = "Kxh6"

    # 3. Black's rook delivers the final blow, protected by the pawn on g6.
    black_move_2 = "Rh4#"

    print(black_move_1)
    print(white_move_1)
    print(black_move_2)

solve_chess_mate()