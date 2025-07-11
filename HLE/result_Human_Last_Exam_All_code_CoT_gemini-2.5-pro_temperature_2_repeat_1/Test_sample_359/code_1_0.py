# Note: This script uses the 'python-chess' library.
# If you don't have it installed, you can install it by running: pip install python-chess

import chess

def solve_koth_puzzle():
    """
    Analyzes the given King of the Hill chess FEN to find how quickly White can win.
    """
    # The FEN string describes the board state.
    fen = "8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43"

    # We use the python-chess library to create a board object from the FEN.
    board = chess.Board(fen)

    # These are the four central squares that lead to a King of the Hill victory.
    koth_winning_squares = [chess.D4, chess.E4, chess.D5, chess.E5]
    koth_square_names = {chess.square_name(s) for s in koth_winning_squares}

    print("--- King of the Hill Puzzle Analysis ---")
    print(f"FEN: {fen}")
    print(f"The central winning squares are: {sorted(list(koth_square_names))}.")

    # Assuming optimal play, White will look for the fastest possible win.
    # We will first check if a win in 1 move is possible.
    print("\nChecking for a win in 1 move for White...")

    winning_move_found = False
    
    # We look at all of White's legal moves.
    for move in board.legal_moves:
        # We are interested in moves made by the King.
        if board.piece_type_at(move.from_square) == chess.KING:
            # Check if the King's destination is a winning central square.
            if move.to_square in koth_winning_squares:
                from_square = chess.square_name(move.from_square)
                to_square = chess.square_name(move.to_square)
                # board.san(move) provides the move in Standard Algebraic Notation.
                move_notation = board.san(move)
                
                print(f"\nConclusion: A winning move has been found!")
                print(f"White can play the move '{move_notation}'.")
                print(f"This moves the King from {from_square} to the central square {to_square}.")
                print("This fulfills the King of the Hill victory condition, so White wins immediately.")
                
                # The number of moves for White to win is 1.
                moves_to_win = 1
                
                # The problem asks to output the number(s) in the final equation.
                # Here, the equation is simple: the number of moves is 1.
                print("\n--- Final Answer ---")
                print(f"Number of moves for White to win = {moves_to_win}")
                
                winning_move_found = True
                break  # We found the fastest win, so we can stop searching.

    if not winning_move_found:
        print("\nNo win in 1 move was found for White. A deeper search would be needed.")

# Execute the analysis function.
solve_koth_puzzle()