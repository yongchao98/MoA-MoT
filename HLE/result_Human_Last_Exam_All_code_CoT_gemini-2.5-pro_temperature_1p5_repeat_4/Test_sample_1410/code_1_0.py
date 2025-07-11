def solve_chess_puzzle():
    """
    This function explains and solves the chess puzzle presented by the FEN:
    4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20
    
    The script will print the move-by-move solution for the forced checkmate.
    """
    
    print("The solution involves a forced checkmate for White.")
    print("Here is the winning sequence, which constitutes the 'equation' of moves:")
    
    # The sequence of moves leading to mate.
    moves = [
        ("Ng6+", "hxg6"),
        ("Qxg6+", "Kh8"),
        ("Qxh6#", "")
    ]
    
    num_white_moves = len(moves)
    
    # Printing the "equation" as requested, with each move number.
    for i, move_pair in enumerate(moves):
        move_number = i + 1
        white_move, black_move = move_pair
        if black_move:
            print(f"{move_number}. {white_move} {black_move}")
        else:
            # The last move is checkmate, no black move follows.
            print(f"{move_number}. {white_move}")
            
    print("\nThis sequence is forced. Let's break down why:")
    print("1. Ng6+: White checks the Black king. The Black pawn on f7 cannot capture because it's pinned by the White bishop on h5. The king cannot move to h8 for the same reason. Black's only legal move is hxg6.")
    print("2. Qxg6+: White checks again. Black's only legal move is to move the king to h8.")
    print("3. Qxh6#: White delivers checkmate. The king is attacked and has no legal moves.")
    
    print("\nThus, assuming optimal play, the number of moves in which White can win is:")
    print(num_white_moves)

solve_chess_puzzle()
<<<3>>>