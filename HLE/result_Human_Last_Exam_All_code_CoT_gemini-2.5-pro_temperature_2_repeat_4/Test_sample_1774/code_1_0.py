import chess
import chess.engine
import sys
import os

def find_stockfish():
    """Tries to find a UCI-compatible chess engine executable."""
    # List of common names and paths for the Stockfish engine
    stockfish_names = [
        "stockfish",
        "stockfish.exe",
        "/usr/games/stockfish",
        "/usr/local/bin/stockfish",
        "/opt/homebrew/bin/stockfish",
    ]
    
    # First, check if stockfish is in the system's PATH
    for path in os.environ.get("PATH", "").split(os.pathsep):
        for name in ["stockfish", "stockfish.exe"]:
            full_path = os.path.join(path, name)
            if os.path.exists(full_path) and os.access(full_path, os.X_OK):
                return full_path

    # If not in PATH, check the other common hardcoded paths
    for path in stockfish_names:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path
            
    return None

def solve_chess_puzzle():
    """
    Sets up the chess board, analyzes it with an engine, and prints the mate sequence.
    """
    try:
        import chess
        import chess.engine
    except ImportError:
        print("Error: The 'python-chess' library is required but not installed.")
        print("Please install it by running: pip install python-chess")
        return

    engine_path = find_stockfish()
    if not engine_path:
        print("---")
        print("Error: Stockfish chess engine not found.")
        print("Please install Stockfish and ensure it's in your system's PATH.")
        print("Installation commands:")
        print("  - On Debian/Ubuntu: sudo apt-get update && sudo apt-get install stockfish")
        print("  - On macOS (with Homebrew): brew install stockfish")
        print("---")
        return

    try:
        engine = chess.engine.SimpleEngine.popen_uci(engine_path)
    except (FileNotFoundError, PermissionError) as e:
        print(f"Error: Could not start the chess engine at '{engine_path}'.")
        print(f"Reason: {e}")
        return

    # FEN string representing the given board state.
    # Board setup: White to move.
    fen = "rn3rk1/pbppq1pp/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1"
    board = chess.Board(fen)

    print("Analyzing the position to find the forced mate...")

    try:
        # Analyze the board for a specific time or depth to find the mate.
        # Depth is generally better for finding forced mates.
        info = engine.analyse(board, chess.engine.Limit(depth=20))
        score = info.get("score").white()

        if score.is_mate():
            mate_in_x = score.moves
            pv = info.get("pv")

            # Create a board copy to play out the moves for display.
            display_board = board.copy()
            move_list_str = []
            
            for i, move in enumerate(pv):
                if display_board.turn == chess.WHITE:
                    move_number = display_board.fullmove_number
                    san_move = display_board.san(move)
                    move_list_str.append(f"{move_number}. {san_move}")
                else:
                    san_move = display_board.san(move)
                    move_list_str[-1] += f" {san_move}"
                
                display_board.push(move)
            
            # Add the checkmate symbol to the final move
            if display_board.is_checkmate():
                 if not move_list_str[-1].endswith('#'):
                    move_list_str[-1] += "#"

            print(f"\nSolution found: White can deliver checkmate in {mate_in_x} moves.")
            print("The winning move sequence is:")
            print(" ".join(move_list_str))
            
            print(f"\n\n<<< {mate_in_x} >>>")

        else:
            print("No forced mate sequence was found within the analysis depth.")
            print(f"The current evaluation is: {score}")

    except Exception as e:
        print(f"An error occurred during analysis: {e}")

    finally:
        engine.quit()

# Execute the solver function
solve_chess_puzzle()