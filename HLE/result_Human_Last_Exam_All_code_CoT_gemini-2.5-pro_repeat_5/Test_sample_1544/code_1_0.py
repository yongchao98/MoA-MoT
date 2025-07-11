import chess
import chess.engine
import sys
import os

def find_best_chess_move():
    """
    Analyzes a chess position using the Stockfish engine to find the best move.

    This script requires the Stockfish chess engine. Please make sure it is
    installed on your system and accessible via your system's PATH, or
    specify the path to the executable directly in the 'stockfish_path' variable.
    """
    fen_position = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    
    # --- Find Stockfish ---
    # List of potential paths for Stockfish executable
    # On Windows, you might have "stockfish.exe" in the same folder as the script.
    # On Linux/macOS, it's often in "/usr/games/stockfish" or "/usr/bin/stockfish".
    potential_paths = [
        # For Windows
        "stockfish.exe", 
        # For Linux/macOS
        "/usr/games/stockfish",
        "/usr/bin/stockfish",
        "/usr/local/bin/stockfish"
    ]
    
    stockfish_path = None
    for path in potential_paths:
        if os.path.exists(path):
            stockfish_path = path
            break

    if not stockfish_path:
        print("Error: Stockfish engine not found.")
        print("Please install Stockfish and ensure its executable is in a common path,")
        print("or add the full path to the 'potential_paths' list in this script.")
        return

    # --- Use Stockfish to find the best move ---
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
        
        board = chess.Board(fen_position)
        
        # Let the engine think for a short time (e.g., 1 second)
        result = engine.play(board, chess.engine.Limit(time=1.0))
        
        # Convert the move to Standard Algebraic Notation (e.g., "Qxc5+")
        best_move_san = board.san(result.move)
        
        print(f"The FEN of the position is: {board.fen()}")
        print(f"The best move for White is: {best_move_san}")

        engine.quit()

    except Exception as e:
        print(f"An error occurred while running Stockfish: {e}")
        # Make sure to close the engine if it was opened
        if 'engine' in locals() and engine.is_alive():
            engine.quit()

if __name__ == "__main__":
    find_best_chess_move()