import chess
import chess.engine
import sys

def find_best_move():
    """
    Analyzes a chess position using the Stockfish engine to find the best move.

    This function requires the python-chess library and a local installation of
    the Stockfish chess engine.

    - Install python-chess: pip install python-chess
    - Download Stockfish: https://stockfishchess.org/download/
    """
    # --- Configuration ---
    # Please update this path to point to your Stockfish executable file.
    # Common paths:
    # Linux: "/usr/games/stockfish" or "/usr/bin/stockfish"
    # Windows: r"C:\path\to\stockfish\stockfish.exe"
    # macOS: "/usr/local/bin/stockfish"
    stockfish_path = "/usr/games/stockfish"
    
    fen_position = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    
    # Check for the engine executable
    try:
        # Give a bit more time for deeper analysis in complex endgames.
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
    except FileNotFoundError:
        print(f"Error: Stockfish engine not found at '{stockfish_path}'.")
        print("Please download Stockfish and update the 'stockfish_path' variable in this script.")
        # Exit if the engine can't be found
        if 'pytest' not in sys.modules:
            sys.exit(1)
        else:
            return "Engine not found"


    board = chess.Board(fen_position)
    
    print("Position to analyze:")
    print(board)
    print(f"\nFEN: {fen_position}")
    print("\nThinking...")

    try:
        # Analyze the position for 3 seconds to find the optimal move.
        result = engine.play(board, chess.engine.Limit(time=3.0))
        
        # Convert the move to Standard Algebraic Notation (SAN).
        best_move_san = board.san(result.move)
        
        print("\n--- Result ---")
        print(f"The best move for White is: {best_move_san}")
        
    except Exception as e:
        print(f"An error occurred during chess analysis: {e}")
    finally:
        # Ensure the engine process is terminated.
        engine.quit()

if __name__ == "__main__":
    find_best_move()
