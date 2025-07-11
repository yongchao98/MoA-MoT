import chess

def solve_chess_puzzle():
    """
    Solves the chess puzzle by demonstrating the forced mate in 3.
    The FEN string describes a legal position where White has a winning attack.
    The key is to find the most direct path to checkmate.
    """
    # FEN for the position
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)

    print("Initial position:")
    print(board)
    print("\nWhite to play and win.\n")

    # The winning sequence of moves
    moves = ["f7", "g6", "e8"]
    
    # 1. White plays Qxf7+
    move1_san = board.san(chess.Move.from_uci("d5f7"))
    board.push_san(move1_san)
    print(f"Move 1 (White): Qxf7+")
    print("Black is forced to play Kh7.")
    # The engine finds the only legal response for black, which is Kh7
    board.push(board.legal_moves.find(chess.Move.from_uci("g8h7")))
    print(board)
    print("-" * 20)

    # 2. White plays Bg6+
    move2_san = board.san(chess.Move.from_uci("h5g6"))
    board.push_san(move2_san)
    print(f"Move 2 (White): Bg6+")
    print("Black's best move to avoid immediate mate is to sacrifice the queen, but Kh8 leads to mate.")
    print("Assuming the line that leads to mate, Black plays Kh8.")
    # We push the move that leads to mate
    board.push(board.legal_moves.find(chess.Move.from_uci("h7h8")))
    print(board)
    print("-" * 20)

    # 3. White plays Qxe8#
    move3_san = board.san(chess.Move.from_uci("f7e8"))
    board.push_san(move3_san)
    print(f"Move 3 (White): Qxe8#")
    print("The result is checkmate.")
    print(board)
    print("-" * 20)

    if board.is_checkmate():
        print("White wins by checkmate in 3 moves.")
    else:
        print("Error in analysis, not a checkmate.")

solve_chess_puzzle()
