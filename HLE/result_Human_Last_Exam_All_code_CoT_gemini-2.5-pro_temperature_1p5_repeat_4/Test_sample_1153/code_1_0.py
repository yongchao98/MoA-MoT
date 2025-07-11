import chess
import sys
import os

def find_stockfish_path():
    """Helper function to find the stockfish executable on common paths."""
    # The 'shutil.which' function is the most reliable way to find an executable.
    try:
        from shutil import which
        for exe_name in ['stockfish', 'stockfish.exe', 'stockfish_15_x64_avx2.exe']:
            path = which(exe_name)
            if path:
                return path
    except ImportError:
        # Fallback for older systems or environments without shutil
        for path in ['/usr/games/stockfish', 'stockfish']:
             if os.path.exists(path) and os.access(path, os.X_OK):
                 return path
    return None

def solve_chess_puzzle():
    """
    Analyzes a chess position with Stockfish to find the best move among a set of options.
    This script requires the 'stockfish' and 'python-chess' libraries.
    Install them with: pip install stockfish python-chess
    It also requires the Stockfish engine executable to be available in the system's PATH.
    """
    try:
        from stockfish import Stockfish
    except ImportError:
        print("Required libraries not found. Please run: pip install stockfish python-chess", file=sys.stderr)
        return

    stockfish_path = find_stockfish_path()
    if not stockfish_path:
        print("Stockfish engine not found.", file=sys.stderr)
        print("Please install the Stockfish engine and ensure its location is in your system's PATH.", file=sys.stderr)
        return
        
    try:
        # Initialize Stockfish with a reasonable depth for a puzzle
        stockfish = Stockfish(path=stockfish_path, depth=18, parameters={"Threads": 2, "Minimum Thinking Time": 1200})
    except Exception as e:
        print(f"Failed to initialize Stockfish: {e}", file=sys.stderr)
        return

    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    # Dictionary of SAN moves for each answer choice
    moves_san = {
        "A": "a8=Q",
        "B": "Nc5",
        "C": "Kd4",
        "D": "Kf2",
        "E": "Nf4",
        "F": "b4",
    }
    
    print(f"Analyzing FEN: {fen}\n")
    
    results = {}

    for key, san_move in moves_san.items():
        # Create a fresh board state for each move analysis
        temp_board = board.copy()
        
        try:
            # push_san finds the move from SAN and applies it. It raises an exception for illegal moves.
            move = temp_board.push_san(san_move)
            
            # Set the engine to the initial position and make the move
            stockfish.set_fen_position(fen)
            stockfish.make_moves_from_current_position([move.uci()])
            
            # The evaluation is now from Black's perspective
            eval_data = stockfish.get_evaluation()
            
            eval_type = eval_data["type"]
            eval_value = eval_data["value"] # This value is for Black

            # Convert to a clear, readable string from White's perspective
            if eval_type == "mate":
                # A negative mate value like -7 means Black is mated in 7, which is a win for White.
                score_str = f"White has a mate in {abs(eval_value)}."
            else: # cp (centipawns)
                # A negative cp value for Black means White has an advantage. We flip the sign.
                white_advantage_cp = -eval_value
                pawn_advantage = white_advantage_cp / 100.0
                score_str = f"White's advantage is {pawn_advantage:+.2f} pawns."
            
            print(f"Move {key} ({san_move}): {score_str}")

        except ValueError:
            # This handles cases where board.push_san fails because the move is illegal
            print(f"Move {key} ({san_move}): This move is illegal.")
        except Exception as e:
            print(f"Move {key} ({san_move}): Could not be evaluated. Error: {e}")

solve_chess_puzzle()