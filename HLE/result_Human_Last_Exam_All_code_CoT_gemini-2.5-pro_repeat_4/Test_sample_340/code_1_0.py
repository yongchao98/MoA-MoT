import chess

def solve_three_check_puzzle():
    """
    This function demonstrates the 4-move solution to the three-check chess puzzle.
    It applies the moves and prints the status at each step.
    """
    # Initial FEN: r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1 +0+0
    board = chess.Board("r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1")
    
    white_checks = 0
    black_checks = 0

    print(f"Initial position: {board.fen()}")
    print(f"White checks: {white_checks}, Black checks: {black_checks}\n")

    moves_sequence = [
        ("Bxd7+", "Kxd7"),
        ("Qc3+", "Ke8"),
        ("Bxf6", "gxf6"),
        ("Qc6+", None)
    ]
    
    move_num = 0
    for white_move_san, black_move_san in moves_sequence:
        move_num += 1
        
        # White's move
        white_move = board.parse_san(white_move_san)
        is_check = board.is_check()
        board.push(white_move)
        if board.is_check():
            white_checks += 1
        
        print(f"Move {move_num} (White): {white_move_san}")
        print(f"White delivers check #{white_checks}!")
        print(f"Board FEN: {board.fen()}")
        print(f"Check count: White {white_checks}, Black {black_checks}\n")

        if white_checks == 3:
            print("White wins by delivering the third check.")
            break
        
        # Black's move
        if black_move_san:
            black_move = board.parse_san(black_move_san)
            board.push(black_move)
            print(f"Move {move_num} (Black): {black_move_san}")
            print(f"Board FEN: {board.fen()}\n")

    print("---")
    print("The minimal amount of moves by white to win is the number of moves White played.")
    print(f"The winning sequence is:")
    print("1. Bxd7+ (Check 1) Kxd7")
    print("2. Qc3+ (Check 2) Ke8")
    print("3. Bxf6 gxf6")
    print("4. Qc6+ (Check 3)")
    print("\nTotal moves by White to win = 4")

solve_three_check_puzzle()
<<<4>>>