import chess
import chess.engine
import sys
import os

def find_stockfish_executable():
    """
    Tries to find a Stockfish executable in common locations or from an
    environment variable.
    """
    # Check for STOCKFISH_PATH environment variable
    stockfish_path = os.getenv("STOCKFISH_PATH")
    if stockfish_path and os.path.isfile(stockfish_path):
        return stockfish_path

    # Common paths for different operating systems
    paths = [
        # Linux
        "/usr/games/stockfish",
        "/usr/bin/stockfish",
        # macOS (via Homebrew)
        "/usr/local/bin/stockfish",
        "/opt/homebrew/bin/stockfish",
        # Windows
        "C:\\Program Files\\Stockfish\\stockfish.exe",
        "stockfish.exe",
        # In current directory
        "./stockfish"
    ]

    for path in paths:
        if os.path.isfile(path):
            return path
            
    return None

def solve_mate_in_x(fen_string: str):
    """
    Analyzes a FEN string to find a forced mate and returns the number of moves.
    """
    stockfish_executable = find_stockfish_executable()
    if not stockfish_executable:
        print("Error: Stockfish engine not found.", file=sys.stderr)
        print("Please install Stockfish and add it to your system's PATH,", file=sys.stderr)
        print("or set the STOCKFISH_PATH environment variable.", file=sys.stderr)
        sys.exit(1)

    engine = None
    try:
        # Start the chess engine
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_executable)
        
        # Load the board position from the FEN string
        board = chess.Board(fen_string)
        
        # Analyze the position for a fixed time to find the best line
        # A 15-second analysis is more than enough for most tactical puzzles.
        info = engine.analyse(board, chess.engine.Limit(time=15.0))
        
        # The score is from the perspective of the side to move (White in this case)
        score = info.get("score")
        
        if score is None:
            print("Could not retrieve a score from the engine.", file=sys.stderr)
            return

        white_score = score.white()
        
        # Check if the score indicates a forced checkmate
        if white_score.is_mate():
            moves_to_mate = white_score.moves
            # Print the number of moves to mate as requested
            print(moves_to_mate)
        else:
            # If no mate is found, print a message indicating so.
            print("No forced checkmate found for White in the given position.")

    except (FileNotFoundError, PermissionError) as e:
        print(f"Error accessing Stockfish: {e}", file=sys.stderr)
        print("Please ensure the Stockfish executable path is correct and has execute permissions.", file=sys.stderr)
    except Exception as e:
        print(f"An unexpected error occurred: {e}", file=sys.stderr)
    finally:
        # Ensure the engine is always closed properly
        if engine:
            engine.quit()

if __name__ == "__main__":
    # The FEN from the problem description
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    solve_mate_in_x(fen)