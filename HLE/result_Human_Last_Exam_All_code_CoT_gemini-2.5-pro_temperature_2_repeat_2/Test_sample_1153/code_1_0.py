import chess
import chess.engine
import sys
import os

def find_stockfish():
    """Finds the Stockfish executable in common locations."""
    # List of common paths for the Stockfish executable
    paths = [
        "/usr/games/stockfish",
        "/usr/bin/stockfish",
        "/opt/homebrew/bin/stockfish",
        "/usr/local/bin/stockfish",
        "stockfish",
        "./stockfish",
    ]
    # For Windows, also check for .exe and common installation directories
    if sys.platform == "win32":
        paths.extend([
            "stockfish.exe",
            "./stockfish.exe",
            "C:\\Program Files\\Stockfish\\stockfish.exe",
        ])

    for path in paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            return path
    return None

def solve_chess_puzzle():
    """
    Analyzes a chess position with Stockfish to find the best move among a list of choices.
    """
    stockfish_path = find_stockfish()
    if not stockfish_path:
        print("---")
        print("Error: Stockfish engine not found.")
        print("Please install Stockfish and add it to your system's PATH,")
        print("or place the executable in the same folder as this script.")
        print("---")
        return

    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
    except (FileNotFoundError, PermissionError):
        print(f"Error: Could not launch stockfish at: {stockfish_path}")
        return

    # FEN for the puzzle
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    # Candidate moves from the question
    moves = {
        "A": "a8=Q",
        "B": "Nc5",
        "C": "Kd4",
        "D": "Kf2",
        "E": "Nf4",
        "F": "b4",
    }

    results = {}
    print(f"Analyzing position: {fen}")
    print("This may take a moment...")
    print("-" * 30)

    # Evaluate each move
    for option, move_san in moves.items():
        try:
            move = board.parse_san(move_san)
            temp_board = board.copy()
            temp_board.push(move)
            
            # Analyze the position for a short duration
            # Increase time for more accurate results
            info = engine.analyse(temp_board, chess.engine.Limit(time=1.0))
            
            # Get the score from White's perspective
            score = info["score"].white()

            # Format the score for printing
            if score.is_mate():
                eval_str = f"Mate in {score.mate()}"
                # A winning mate is given a very high value for sorting
                eval_val = 10000 - score.mate()
            else:
                # Convert centipawns to a standard pawn advantage score
                pawn_advantage = score.score() / 100.0
                eval_str = f"+{pawn_advantage:.2f}"
                eval_val = pawn_advantage
            
            results[option] = (move_san, eval_str, eval_val)
        
        except (chess.InvalidMoveError, chess.AmbiguousMoveError) as e:
            results[option] = (move_san, f"Error: {e}", -999)
    
    engine.quit()

    if not results:
        print("Could not analyze any moves.")
        return
        
    # Sort results by the original option letter for consistent display
    sorted_options = sorted(results.items())

    print("Evaluation of each candidate move:")
    # This loop fulfills the instruction to output each number in the final equation.
    for option, (move_san, eval_str, _) in sorted_options:
        # We format the output to look like a clear equation
        print(f"({option}) {move_san:<5} evaluation = {eval_str} pawns")

    # Find the best move by comparing the numerical evaluation values
    best_option_key = max(results, key=lambda k: results[k][2])
    best_move_san = results[best_option_key][0]

    print("-" * 30)
    print(f"The best move is ({best_option_key}) {best_move_san}, as it leads to the greatest advantage for White.")


if __name__ == '__main__':
    solve_chess_puzzle()
<<<B>>>