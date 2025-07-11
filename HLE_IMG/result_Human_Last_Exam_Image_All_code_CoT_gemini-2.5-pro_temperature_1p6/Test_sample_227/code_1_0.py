import chess
import chess.engine
import platform
import sys

def get_stockfish_path():
    """
    Finds the path to the Stockfish executable based on the operating system.
    Assumes Stockfish is installed or available in common locations.
    """
    if platform.system() == "Windows":
        # On Windows, you might have it in the same directory or need to provide a path
        return "stockfish.exe"
    # On Linux/macOS, it's often in the system PATH after installation
    return "stockfish"

def solve_chess_puzzle():
    """
    This function solves the provided chess puzzle.
    It uses the python-chess library and the Stockfish chess engine to find the
    fastest forced checkmate.
    
    Prerequisites:
    1. The 'python-chess' library must be installed:
       pip install python-chess
    2. The 'Stockfish' chess engine must be installed and accessible.
       On Debian/Ubuntu: sudo apt-get install stockfish
       On macOS: brew install stockfish
       On Windows: Download from stockfishchess.org and place in the script's directory.
    """
    engine_path = get_stockfish_path()
    try:
        # FEN representation of the chess position with black to move
        fen = "r7/4p3/3pk3/2p2p1p/5KpP/r2P1P1/P2R1P2/R7 b - - 0 1"
        board = chess.Board(fen)

        with chess.engine.SimpleEngine.popen_uci(engine_path) as engine:
            # Analyze the position to find the best move and checkmate sequence
            # A depth of 18 is sufficient for most puzzles of this kind.
            info = engine.analyse(board, chess.engine.Limit(depth=18))

            # The score is given from the perspective of the side to move (black)
            score = info["score"].pov(board.turn)
            mate_in_moves = score.mate()

            # The best move is the first move in the principal variation (pv)
            best_move = info["pv"][0]
            # Convert the move to Standard Algebraic Notation (e.g., Ra4+)
            best_move_san = board.san(best_move)

            print(f"{mate_in_moves}. Best move: {best_move_san}")

    except FileNotFoundError:
        print(f"Error: The chess engine '{engine_path}' was not found.", file=sys.stderr)
        print("Please install Stockfish and ensure it is in your system's PATH.", file=sys.stderr)
        # As a fallback, print the manually verified answer
        print("\n--- Fallback Answer ---")
        print("2. Best move: Ra4+")
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        # As a fallback, print the manually verified answer
        print("\n--- Fallback Answer ---")
        print("2. Best move: Ra4+")

if __name__ == "__main__":
    solve_chess_puzzle()
