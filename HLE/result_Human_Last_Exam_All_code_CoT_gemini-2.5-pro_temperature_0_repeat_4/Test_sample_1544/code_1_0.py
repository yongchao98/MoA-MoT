import chess
import chess.engine
import shutil

def solve_chess_puzzle():
    """
    This script finds the best move for White in a given chess position
    using the Stockfish chess engine.
    """
    # --- Prerequisites ---
    # This script requires the 'python-chess' library and the 'Stockfish' chess engine.
    #
    # To install the library, run this command in your terminal:
    # pip install python-chess
    #
    # To install Stockfish on Debian/Ubuntu/Google Colab, run:
    # sudo apt-get update && sudo apt-get install stockfish
    #
    # On other systems (like Windows or macOS), please download it from
    # https://stockfishchess.org/download/ and ensure the executable is
    # in your system's PATH, or update the 'stockfish_path' variable below.

    # --- Find Stockfish Engine ---
    # The script will try to find the Stockfish executable automatically.
    stockfish_path = shutil.which("stockfish")
    if not stockfish_path:
        # Common path on Linux if not in PATH
        stockfish_path = "/usr/games/stockfish"
        if not shutil.os.path.exists(stockfish_path):
            print("---")
            print("Error: Stockfish engine not found.")
            print("Please install Stockfish and ensure it's in your system's PATH,")
            print("or update the 'stockfish_path' variable in this script with the correct location.")
            print("---")
            return

    engine = None
    try:
        # --- Engine and Board Setup ---
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
        
        # FEN for the given chess position
        fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
        board = chess.Board(fen)

        # --- Analysis ---
        # Ask the engine to analyze the position for 2 seconds.
        # This is usually enough for a clear tactical solution.
        print("Analyzing position... Please wait.")
        result = engine.play(board, chess.engine.Limit(time=2.0))
        best_move = result.move
        
        # Convert the move from UCI format (e.g., f2b2) to SAN format (e.g., Qb2+)
        san_move = board.san(best_move)

        # --- Output ---
        print("\n--- Analysis Complete ---")
        print(f"Position FEN: {board.fen()}")
        print("Assuming optimal play from both sides...")
        print(f"White's best move is: {san_move}")
        print("-------------------------")

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        # --- Cleanup ---
        # Ensure the engine process is terminated
        if engine:
            engine.quit()

# Run the solver
if __name__ == "__main__":
    solve_chess_puzzle()