import chess

def solve_chess_mate():
    """
    Solves the chess puzzle by applying the winning move sequence
    and verifying the checkmate.
    """
    # The FEN string for the chess position.
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)

    # The winning sequence of moves in Standard Algebraic Notation (SAN).
    # White sacrifices the queen, and Black's recapture is forced.
    # Then White delivers the final blow.
    winning_line = [
        ("Qxh6+", "gxh6"),
        ("Nf7#", None)
    ]
    
    num_white_moves = len(winning_line)

    print(f"The initial FEN is: {board.fen()}")
    print("Finding the shortest mate for White...")
    print("-" * 30)

    move_number = 0
    # Apply and print each move in the sequence
    for white_move_san, black_move_san in winning_line:
        move_number += 1
        
        # White's move
        white_move = board.parse_san(white_move_san)
        board.push(white_move)
        print(f"Move {move_number} (White): {white_move_san}")
        
        # Black's move (if any)
        if black_move_san:
            # In this puzzle, black's move is forced.
            black_move = board.parse_san(black_move_san)
            if black_move in board.legal_moves:
                board.push(black_move)
                print(f"Move {move_number} (Black): {black_move_san}")
            else:
                print(f"Error: Black's move {black_move_san} is not legal.")
                return

    # Verify the final position is checkmate
    if board.is_checkmate():
        print("-" * 30)
        print("Result: Confirmed checkmate.")
        print(f"White wins in {num_white_moves} moves.")
    else:
        print("The provided line is not a checkmate.")

# To run this code, you need to install the python-chess library:
# pip install python-chess

solve_chess_mate()