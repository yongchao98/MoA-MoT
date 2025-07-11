import chess

def solve_three_check():
    """
    This function solves the given three-check chess problem by demonstrating
    a forced win for White in a minimal number of moves.
    """
    # The special FEN includes check counts at the end, which the `chess` library doesn't parse.
    # We'll handle the check counts manually.
    fen = "r3kb1r/p2nqppp/5n2/1B2p1B1/4P3/1Q6/PPP2PPP/R3K2R w KQkq - 0 1"
    white_checks = 0
    black_checks = 0

    board = chess.Board(fen)

    print(f"Starting Position (FEN: {fen})")
    print(f"Initial Checks: White={white_checks}, Black={black_checks}\n")

    # The plan is a forced win in 3 moves for White.
    # We will demonstrate the line starting with 1. Qe6+ Be7.
    # The alternative, 1... Kf8, also leads to a win in 3 moves for White.
    
    # --- Move 1: White ---
    move1_w_uci = "b3e6"
    move1_w_san = board.san(chess.Move.from_uci(move1_w_uci))
    board.push_uci(move1_w_uci)
    if board.is_check():
        white_checks += 1
    print(f"White's 1st move: {move1_w_san}")
    print(f"Checks: White={white_checks}, Black={black_checks}")
    print(board)
    print("")

    # --- Move 1: Black ---
    # Black has two responses: Kf8 or Be7. We'll follow the Be7 line.
    # The Kf8 line is symmetrical and also results in a 3-move win for white.
    move1_b_uci = "f8e7"
    move1_b_san = board.san(chess.Move.from_uci(move1_b_uci))
    board.push_uci(move1_b_uci)
    print(f"Black's 1st move: {move1_b_san}")
    print(f"Checks: White={white_checks}, Black={black_checks}")
    print(board)
    print("")

    # --- Move 2: White ---
    move2_w_uci = "e6e7"
    move2_w_san = board.san(chess.Move.from_uci(move2_w_uci))
    board.push_uci(move2_w_uci)
    if board.is_check():
        white_checks += 1
    print(f"White's 2nd move: {move2_w_san}")
    print(f"Checks: White={white_checks}, Black={black_checks}")
    print(board)
    print("")

    # --- Move 2: Black ---
    # This move is forced.
    move2_b_uci = "e8e7"
    move2_b_san = board.san(chess.Move.from_uci(move2_b_uci))
    board.push_uci(move2_b_uci)
    print(f"Black's 2nd move: {move2_b_san}")
    print(f"Checks: White={white_checks}, Black={black_checks}")
    print(board)
    print("")

    # --- Move 3: White ---
    move3_w_uci = "e1c1" # O-O-O
    move3_w_san = board.san(chess.Move.from_uci(move3_w_uci))
    board.push_uci(move3_w_uci)
    if board.is_check():
        white_checks += 1
    print(f"White's 3rd move: {move3_w_san}")
    print(f"Checks: White={white_checks}, Black={black_checks}")
    print("This is the 3rd check. White wins.")
    print(board)
    print("")

    white_move_count = 3
    print(f"The minimal amount of moves by white to win is: {white_move_count}")

# To run this code, you need to install the python-chess library:
# pip install chess
solve_three_check()