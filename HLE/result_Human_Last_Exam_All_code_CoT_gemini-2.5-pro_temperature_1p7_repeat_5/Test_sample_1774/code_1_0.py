import chess
import chess.engine
import sys
import os

def find_stockfish_path():
    """Helper function to find the Stockfish executable."""
    # Common paths for different operating systems
    paths = [
        "stockfish",  # In system PATH
        "/usr/games/stockfish",  # Common on Debian/Ubuntu
        "/usr/bin/stockfish",    # Common on other Linux distros
        "/opt/homebrew/bin/stockfish",  # For Apple Silicon Macs with Homebrew
        "./stockfish.exe" # For Windows in the same directory
    ]
    for path in paths:
        # On Windows, we need to check for .exe
        if sys.platform == "win32" and not path.endswith(".exe"):
            path_exe = path + ".exe"
            if os.path.exists(path_exe):
                return path_exe
        
        # Check if the path is a valid file
        if os.path.exists(path):
             # Check if it's executable
             if os.access(path, os.X_OK):
                 return path
    
    # Try an absolute path for Windows as a last resort
    if sys.platform == "win32":
       if os.path.exists("C:/Stockfish/stockfish.exe"):
          return "C:/Stockfish/stockfish.exe"
          
    return None

def solve_chess_mate():
    """
    Analyzes the given chess position to find the shortest mate for White.
    """
    # 1. Set up the board from the FEN string derived from the problem description.
    fen = "rn3r1k/pbppq1pp/1p2p1b1/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w - - 0 1"
    board = chess.Board(fen)

    # 2. Initialize the Stockfish chess engine.
    stockfish_path = find_stockfish_path()
    if not stockfish_path:
        print("Error: Stockfish engine not found.")
        print("Please install Stockfish and ensure it's in your system PATH,")
        print("or place the executable in the same directory as this script.")
        return

    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
    except Exception as e:
        print(f"Error initializing Stockfish engine: {e}")
        return

    # 3. Analyze the position to find a mate.
    # A shallow search depth is sufficient for such a tactical position.
    try:
        info = engine.analyse(board, chess.engine.Limit(depth=20))
    except Exception as e:
        print(f"An error occurred during engine analysis: {e}")
        engine.quit()
        return

    # 4. Extract and print the results.
    score = info.get("score")
    pv = info.get("pv")  # Principal Variation (best sequence of moves)

    if score and score.is_mate():
        mate_in = score.white().mate()

        if mate_in > 0:
            print(f"White has a forced checkmate in {mate_in} moves.")
            print("\nThe winning sequence is:")
            
            # Build the final string showing the move sequence, like "1. e4 e5 2. Nf3 ..."
            move_log = []
            temp_board = board.copy()
            for i, move in enumerate(pv):
                if temp_board.turn == chess.WHITE:
                    move_log.append(f"{temp_board.fullmove_number}.")
                
                # The python-chess library can generate the move in Standard Algebraic Notation (SAN)
                san_move = temp_board.san(move)
                move_log.append(san_move)
                temp_board.push(move)

            print(" ".join(move_log))
            
            # 5. Output the final answer in the required format.
            print(f"\n<<<{mate_in}>>>")
        else:
            print("The analysis indicates Black is the one with a mate.")
            print("<<<No mate for White>>>")

    else:
        print("No forced mate was found in the analysis.")
        print(f"Current evaluation: {score}")
        print("<<<No mate for White>>>")

    # Clean up the engine process
    engine.quit()

if __name__ == "__main__":
    solve_chess_mate()