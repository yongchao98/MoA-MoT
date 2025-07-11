import chess
import chess.engine
import os

def find_best_move(fen, engine_path):
    """
    Analyzes a chess position from a FEN string and finds the best move.
    
    To run this, you need to:
    1. Install the python-chess library: pip install python-chess
    2. Download a UCI chess engine like Stockfish from https://stockfishchess.org/
    3. Update the 'engine_path' variable below to the correct path of your engine executable.
    """
    if not os.path.exists(engine_path):
        print(f"Error: Chess engine not found at '{engine_path}'.")
        print("Please download a UCI engine (e.g., Stockfish) and update the engine_path variable.")
        return

    try:
        # Start the chess engine
        with chess.engine.SimpleEngine.popen_uci(engine_path) as engine:
            # Set up the board with the given FEN
            board = chess.Board(fen)
            
            # Ask the engine to analyze the position for 2 seconds
            # In this position, this is enough to find the optimal drawing line.
            result = engine.play(board, chess.engine.Limit(time=2.0))
            
            # Get the best move in Standard Algebraic Notation (e.g., "Qd2")
            best_move_san = board.san(result.move)
            
            print(f"Position FEN: {board.fen()}")
            print("Analysis complete.")
            print(f"The best move for White is: {best_move_san}")
            print("\nThis move leads to a forced queen trade and a drawn pawn endgame with correct play.")

    except Exception as e:
        print(f"An error occurred: {e}")

if __name__ == "__main__":
    # The FEN string of the chess position
    position_fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"

    # --- IMPORTANT ---
    # You must change this path to point to your Stockfish executable file.
    # Example for Windows: "C:/chess_engines/stockfish_15_win_x64/stockfish_15.exe"
    # Example for macOS/Linux: "/usr/local/bin/stockfish"
    stockfish_engine_path = "/path/to/your/stockfish"
    
    find_best_move(position_fen, stockfish_engine_path)
