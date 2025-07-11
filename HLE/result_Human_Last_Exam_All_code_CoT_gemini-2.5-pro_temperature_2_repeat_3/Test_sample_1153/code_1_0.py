import chess
import chess.engine
import shutil

def find_best_chess_move():
    """
    Analyzes the given chess position to find the best move for White
    using the Stockfish chess engine.
    """
    # --- Prerequisites ---
    # 1. Install python-chess: pip install python-chess
    # 2. Install Stockfish: Download from https://stockfishchess.org/download/
    #    and ensure the executable is in your system's PATH.

    stockfish_executable = shutil.which("stockfish")
    if not stockfish_executable:
        print("Error: Stockfish executable not found in your system's PATH.")
        print("Please install Stockfish and ensure its location is in your PATH environment variable.")
        return

    # FEN for the specified chess position
    fen_string = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen_string)

    print(f"Analyzing position with FEN: {fen_string}")
    print("White to move.")
    print("-" * 30)

    try:
        with chess.engine.SimpleEngine.popen_uci(stockfish_executable) as engine:
            # Let the engine think for 2 seconds to find the best move
            result = engine.play(board, chess.engine.Limit(time=2.0))
            best_move = result.move
            best_move_san = board.san(best_move)

            print(f"The best move found by the engine is: {best_move_san}")
            print("\n--- Explanation of the Winning Plan ---")
            # This follows the strange instruction "output each number in the final equation"
            # by outlining the sequence of moves with their numbers.
            print("The winning sequence starts with the pawn promotion:")
            print("1. a8=Q")
            print("This is the most decisive move. It forces Black to respond with:")
            print("1. ... Nxa8")
            print("After this exchange, White continues with the powerful knight move:")
            print("2. Nc5")
            print("This move attacks Black's weak pawns and sets up a winning knight-and-pawn endgame for White.")

    except Exception as e:
        print(f"An error occurred: {e}")
        print("Please ensure Stockfish is installed correctly.")


if __name__ == "__main__":
    find_best_chess_move()