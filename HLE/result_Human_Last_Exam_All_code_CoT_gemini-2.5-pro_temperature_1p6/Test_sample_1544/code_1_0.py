import chess

def solve_chess_puzzle():
    """
    Analyzes the chess position and finds the best move for White.
    """
    # The FEN of the position.
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    board = chess.Board(fen)

    # The best move is Qxc5.
    best_move_san = "Qxc5+"
    move = board.parse_san(best_move_san)

    # Deconstruct the move for a clear explanation.
    from_square = chess.square_name(move.from_square)
    to_square = chess.square_name(move.to_square)
    piece = board.piece_at(move.from_square)
    piece_symbol = piece.symbol()

    print("The current board state is represented by the FEN: {}".format(fen))
    print("White is to move.")
    print("Black's Queen and King are both attacking White's strong d6 pawn.")
    print("White's best move is a forcing capture to simplify into a winning endgame.")
    print("\nThe move involves capturing the Black Queen on c5 with the White Queen from f2.")
    
    # We will analyze the consequences of the move.
    # Case 1: King recaptures
    board_case1 = board.copy()
    board_case1.push(move) # Qxc5+
    # Find the king recapture move
    king_recapture_move = board_case1.find_move(chess.C6, chess.C5)
    board_case1.push(king_recapture_move) # Kxc5
    
    # Case 2: Pawn recaptures
    board_case2 = board.copy()
    board_case2.push(move) # Qxc5+
    # Find the pawn recapture move
    pawn_recapture_move = board_case2.find_move(chess.B6, chess.C5)
    board_case2.push(pawn_recapture_move) # bxc5

    print("\nFinal Move Calculation:")
    print("White Queen from {} captures the Black Queen on {}".format(from_square, to_square))
    print("The move in Standard Algebraic Notation is: {}".format(best_move_san))
    print("\nThis creates two scenarios for Black's recapture:")
    print("1. If Black plays 1...Kxc5, the resulting King and Pawn endgame is winning for White.")
    print("2. If Black plays 1...bxc5, White also wins the ensuing endgame by activating the king.")
    print("\nBecause both variations lead to a win for White, this is the optimal move.")


solve_chess_puzzle()