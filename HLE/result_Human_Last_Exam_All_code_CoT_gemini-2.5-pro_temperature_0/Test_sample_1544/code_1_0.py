import chess
import chess.engine
import sys
import os

def find_best_move():
    """
    Analyzes a chess position with Stockfish to find the best move.
    """
    # --- Configuration ---
    # IMPORTANT: Update this path to where you have the Stockfish executable
    # Common paths are checked, but you may need to set it manually.
    # Windows example: "C:/Users/YourUser/Downloads/stockfish/stockfish.exe"
    # macOS/Linux example: "/usr/local/bin/stockfish" or "./stockfish"
    STOCKFISH_PATH = "stockfish" # Default name, assumes it's in your system's PATH

    # --- FEN for the chess position ---
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    board = chess.Board(fen)

    # --- Find and launch the engine ---
    engine = None
    try:
        # Check common paths for Stockfish if the default one fails
        possible_paths = [
            STOCKFISH_PATH,
            "/usr/games/stockfish", # Common on Linux
            "C:\\Program Files\\Stockfish\\stockfish.exe", # Common on Windows
            os.path.join(os.path.expanduser("~"), "Downloads", "stockfish", "stockfish.exe"),
            os.path.join(os.path.expanduser("~"), "Downloads", "stockfish")
        ]
        
        found_path = None
        for path in possible_paths:
            try:
                # Test if the path is a valid executable
                with chess.engine.SimpleEngine.popen_uci(path) as test_engine:
                    found_path = path
                break
            except (FileNotFoundError, PermissionError):
                continue

        if not found_path:
            print("---")
            print("Error: Stockfish engine not found.")
            print("Please download Stockfish from https://stockfishchess.org/download/")
            print("Then, update the STOCKFISH_PATH variable in this script to the correct file path.")
            print("---")
            sys.exit(1)

        engine = chess.engine.SimpleEngine.popen_uci(found_path)

        # --- Analyze the position ---
        # We give the engine 2 seconds to think.
        result = engine.play(board, chess.engine.Limit(time=2.0))
        best_move = result.move

        # --- Print the result ---
        # The move is converted to Standard Algebraic Notation (e.g., "Qxc5+").
        san_move = board.san(best_move)
        print(f"The best move for White is: {san_move}")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure python-chess is installed (`pip install python-chess`) and Stockfish is correctly configured.")

    finally:
        # --- Clean up ---
        if engine:
            engine.quit()

if __name__ == "__main__":
    find_best_move()
<<<Qxc5+>>>