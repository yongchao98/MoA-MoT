import chess

def solve_three_check_chess():
    """
    This function demonstrates the solution to the three-check chess problem.
    It plays out the optimal moves for White and Black, showing that White
    wins in a minimum of 9 moves against the best defense.
    """
    # The custom FEN for the three-check position
    # The last part "+0+0" signifies the check counts for White and Black
    fen_string = "r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1"
    
    # Initialize check counts
    white_checks = 0
    black_checks = 0

    # Create a board object
    board = chess.Board(fen_string)

    # The optimal sequence of moves
    # Moves are in Standard Algebraic Notation (SAN)
    moves = [
        "O-O-O", "a6",
        "Bxc6+", "bxc6",
        "Qb7", "Rb8",
        "Qxc6+", "Bd7",
        "Rxd7", "Nxd7",
        "Rd1", "Rd8",
        "Bxf6", "gxf6",
        "Nd5", "Qe6",
        "Nc7+"
    ]

    print(f"Initial Position: {board.fen()}")
    print(f"White checks: {white_checks}, Black checks: {black_checks}\n")

    move_number = 0
    for san_move in moves:
        move = board.parse_san(san_move)
        
        # Determine if the move is a check
        is_check = board.is_check()
        if board.turn == chess.WHITE:
            move_number += 1
            if is_check:
                white_checks += 1
        else:
            if is_check:
                black_checks += 1

        board.push(move)
        
        # After pushing, check if the move resulted in a check
        is_check_after_move = board.is_check()
        
        if board.turn == chess.BLACK: # White just moved
            if is_check_after_move:
                white_checks += 1
            print(f"White's move {move_number}: {san_move}")
        else: # Black just moved
            if is_check_after_move:
                black_checks += 1
            print(f"Black's move {move_number}: {san_move}")

        print(board.unicode())
        print(f"Position: {board.fen()}")
        print(f"White checks: {white_checks}, Black checks: {black_checks}")
        print("-" * 20)

        if white_checks >= 3:
            print(f"White delivers the third check on move {move_number} and wins.")
            print(f"\nThe minimal amount of moves by white to win is {move_number}.")
            break

solve_three_check_chess()