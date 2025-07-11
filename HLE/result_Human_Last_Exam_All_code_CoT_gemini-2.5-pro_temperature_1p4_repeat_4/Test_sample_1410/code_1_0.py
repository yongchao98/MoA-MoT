import chess
import chess.engine
import sys
import os
import shutil

def find_stockfish():
    """Finds the stockfish executable in common locations or PATH."""
    # Check if stockfish is in the system's PATH
    stockfish_path = shutil.which("stockfish")
    if stockfish_path:
        return stockfish_path

    # Check common installation paths for different OS
    paths_to_check = [
        "/usr/games/stockfish",       # Common on Debian/Ubuntu
        "/usr/bin/stockfish",         # Common on other Linux distros
        "/opt/homebrew/bin/stockfish",# For Homebrew on macOS
    ]
    if sys.platform == "win32":
        # Common paths on Windows
        paths_to_check.extend([
            os.path.expandvars(r"%ProgramFiles%\Stockfish\stockfish.exe"),
            os.path.expandvars(r"%ProgramFiles(x86)%\Stockfish\stockfish.exe"),
            "stockfish.exe" # If it's in the current directory
        ])

    for path in paths_to_check:
        if os.path.exists(path):
            return path
            
    return None

def solve_mate_in_x():
    """
    Analyzes the given FEN to find the number of moves to a forced checkmate.
    """
    FEN = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    
    stockfish_path = find_stockfish()
    if not stockfish_path:
        print("Error: Stockfish engine not found.", file=sys.stderr)
        print("Please install Stockfish and ensure it's in your system's PATH or a standard location.", file=sys.stderr)
        sys.exit(1)

    engine = None
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
        board = chess.Board(FEN)

        # Analyze the position to find a mate. A depth search is effective for this.
        # Depth 20 is sufficient for most short tactical puzzles.
        info = engine.analyse(board, chess.engine.Limit(depth=20))

        # The score is from the perspective of the side to move (White).
        score = info.get("score").white()

        if score.is_mate():
            num_moves_to_win = score.moves
            # Per the instructions, we output the number for the final "equation".
            # For a "mate in X" problem, the solution is the integer X.
            # The simplest equation is Result = X. We print the number X.
            print(num_moves_to_win)
        else:
            # This case shouldn't be reached for the given FEN.
            print("A forced checkmate was not found at the specified analysis depth.", file=sys.stderr)
            
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
    finally:
        if engine:
            engine.quit()

if __name__ == "__main__":
    solve_mate_in_x()