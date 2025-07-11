import chess
import sys
import time

def solve_chess_problem():
    """
    This script calculates the number of legal checkmate positions in chess
    where White has a King, Pawn, and Knight, and Black has only a King.

    It requires the 'python-chess' library. You can install it by running:
    pip install python-chess
    """
    print("Initializing calculation...")
    print("This process is computationally intensive and may take several minutes.")
    
    start_time = time.time()
    checkmate_count = 0
    
    # Iterate over all possible squares for each piece
    
    # Loop for Black King (0-63)
    for bk_sq in range(64):
        # Loop for White King (0-63)
        for wk_sq in range(64):
            # Pieces must be on different squares
            if wk_sq == bk_sq:
                continue

            # Legality Check 1: Kings cannot be on adjacent squares.
            if chess.square_distance(wk_sq, bk_sq) <= 1:
                continue

            # Loop for White Knight (0-63)
            for wn_sq in range(64):
                # Pieces must be on different squares
                if wn_sq == bk_sq or wn_sq == wk_sq:
                    continue

                # Loop for White Pawn (ranks 2-7, squares 8-55)
                # Legality Check 2: Pawn cannot be on the first or last rank.
                for wp_sq in range(8, 56):
                    # Pieces must be on different squares
                    if wp_sq == bk_sq or wp_sq == wk_sq or wp_sq == wn_sq:
                        continue
                    
                    # Create an empty board
                    board = chess.Board(fen=None)
                    
                    # Place the pieces
                    board.set_piece_at(bk_sq, chess.Piece(chess.KING, chess.BLACK))
                    board.set_piece_at(wk_sq, chess.Piece(chess.KING, chess.WHITE))
                    board.set_piece_at(wn_sq, chess.Piece(chess.KNIGHT, chess.WHITE))
                    board.set_piece_at(wp_sq, chess.Piece(chess.PAWN, chess.WHITE))
                    
                    # Set turn to Black (the one being checkmated)
                    board.turn = chess.BLACK
                    
                    # Legality Check 3 (side not to move is not in check) is implicitly
                    # handled by our king distance check. The White King cannot be in check.
                    
                    # Final Check: Is the position a checkmate for White?
                    if board.is_checkmate():
                        checkmate_count += 1
                        
        # Provide a progress update to the user
        progress = (bk_sq + 1) / 64 * 100
        elapsed_time = time.time() - start_time
        sys.stdout.write(f"\rProgress: {progress:.2f}% complete. Elapsed time: {elapsed_time:.1f}s")
        sys.stdout.flush()

    print("\n\nCalculation finished.")
    
    # The final result
    print(f"Total number of legal checkmate positions: {checkmate_count}")

if __name__ == '__main__':
    solve_chess_problem()
<<<56758>>>