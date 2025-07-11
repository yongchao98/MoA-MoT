import chess
import chess.engine
import sys
import os

def find_best_move(fen_string):
    """
    Analyzes a chess position using the Stockfish engine to find the best move.

    This function requires the Stockfish chess engine executable to be available.
    Please install Stockfish from https://stockfishchess.org/download/
    and ensure its location is in your system's PATH or provide the direct path below.
    """
    # --- IMPORTANT ---
    # Update this path to your Stockfish executable if it's not in the system's PATH
    # Examples:
    # Windows: "C:/Users/user/Downloads/stockfish/stockfish.exe"
    # macOS (Homebrew): "/opt/homebrew/bin/stockfish"
    # Linux: "/usr/games/stockfish"
    stockfish_path = "stockfish"

    try:
        # Check if the provided path is a file, otherwise search in PATH
        if not os.path.isfile(stockfish_path):
             which_stockfish = os.popen('which stockfish').read().strip()
             if which_stockfish:
                 stockfish_path = which_stockfish
             else:
                 raise FileNotFoundError

        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
    except FileNotFoundError:
        print("Error: Stockfish engine not found.", file=sys.stderr)
        print("Please install Stockfish and add it to your system's PATH, or update the 'stockfish_path' variable in this script.", file=sys.stderr)
        return

    board = chess.Board(fen_string)

    print("Analyzing the following position for White:")
    print(board)
    print("\nThis may take a few seconds...")

    # Analyze the position for 3 seconds to find the best move.
    try:
        result = engine.analyse(board, chess.engine.Limit(time=3.0))
        # The 'pv' (principal variation) key holds the sequence of best moves.
        best_move = result['pv'][0]
        # The 'score' key holds the evaluation of the position.
        score = result['score'].white()
    except Exception as e:
        print(f"An error occurred during engine analysis: {e}", file=sys.stderr)
        engine.quit()
        return

    # Convert the move to Standard Algebraic Notation for easy reading.
    best_move_san = board.san(best_move)
    engine.quit()

    print("-" * 30)
    print(f"Analysis Complete. Assuming optimal play from both sides:")
    print(f"White's best move is: {best_move_san}")

    # The 'equation' part can be interpreted as the engine's evaluation score.
    # This number quantifies the advantage.
    if score.is_mate():
        # A mate score is the most definitive 'number' in chess analysis.
        print(f"This move leads to a forced checkmate in {score.mate()} moves.")
    else:
        # A centipawn score shows the material advantage. +100 is like being up one pawn.
        pawn_advantage = score.score() / 100.0
        print(f"The engine's evaluation after this move is a decisive advantage of {pawn_advantage} pawns.")

if __name__ == "__main__":
    # The FEN string from your request
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    find_best_move(fen)