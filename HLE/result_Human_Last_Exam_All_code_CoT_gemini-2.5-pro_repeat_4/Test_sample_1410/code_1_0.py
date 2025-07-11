import chess
import chess.engine
import sys
import os

def solve_chess_puzzle():
    """
    Analyzes a chess position to find a forced mate.
    """
    # --- Configuration ---
    # IMPORTANT: Update this path to your Stockfish executable.
    # On Windows, it might be r"C:\Users\YourUser\Downloads\stockfish\stockfish.exe"
    # On Linux/macOS, it might be "/usr/games/stockfish" or a local path like "./stockfish"
    # We will try to find a common path, but you may need to change it.
    
    # List of potential paths for the Stockfish engine
    potential_paths = [
        # For Linux/macOS in common locations
        "/usr/games/stockfish",
        "/usr/local/bin/stockfish",
        # For local directory (e.g., if you downloaded it to the same folder)
        "stockfish",
        "stockfish.exe" 
    ]
    
    STOCKFISH_PATH = None
    for path in potential_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            STOCKFISH_PATH = path
            break

    if STOCKFISH_PATH is None:
        print("---")
        print("Error: Stockfish engine not found or is not executable.")
        print("Please download the Stockfish chess engine, place it in this script's directory,")
        print("or update the `potential_paths` list with the correct path to your stockfish executable.")
        print("---")
        return

    # FEN for the specified chess position
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"

    try:
        engine = chess.engine.SimpleEngine.popen_uci(STOCKFISH_PATH)
    except Exception as e:
        print(f"Failed to start the engine. Error: {e}")
        return

    board = chess.Board(fen)
    print(f"Analyzing chess position: {fen}")
    print("White to move.")
    print("-" * 20)

    try:
        # Analyze for a sufficient depth to find the mate. Depth 20 is plenty for a mate-in-5.
        info = engine.analyse(board, chess.engine.Limit(depth=20))
    except Exception as e:
        print(f"An error occurred during analysis: {e}")
        engine.quit()
        return

    # Get the score from White's perspective
    score = info.get("score").white()

    if score.is_mate():
        mate_in_moves = score.mate()
        print(f"Result: White can force a checkmate in {mate_in_moves} moves.")

        # The phrase "output each number in the final equation" is interpreted as showing
        # the full sequence of moves leading to the checkmate.
        print("\nThe winning equation (sequence of moves) is:")
        
        pv = info.get("pv")
        if pv:
            temp_board = board.copy()
            # Iterate through the principal variation to show the moves
            for i in range(mate_in_moves * 2 - 1):
                if i >= len(pv):
                    break
                move = pv[i]
                if temp_board.turn == chess.WHITE:
                    # Print move number for White's turn
                    print(f"{temp_board.fullmove_number}. {temp_board.san(move)}", end=' ')
                else:
                    # Print Black's move and a newline
                    print(f"{temp_board.san(move)}")
                temp_board.push(move)
            
            # Add a final newline if the sequence ends on a White move
            if temp_board.turn == chess.BLACK:
                print()
            
            if temp_board.is_checkmate():
                print("\nFinal position is checkmate.")

    else:
        print("Result: No forced checkmate found for White in the given depth.")
        print(f"The engine's evaluation is: {score}")

    engine.quit()

if __name__ == "__main__":
    solve_chess_puzzle()