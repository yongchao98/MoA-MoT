import chess
import chess.engine
import os

def find_best_chess_move(fen_string):
    """
    Analyzes a chess position using the Stockfish engine to find the best move.

    This function requires the Stockfish engine executable to be available in the
    system's PATH or in the same directory as the script.
    """
    try:
        # On Windows, you might need to provide the full path to stockfish.exe
        # e.g., chess.engine.SimpleEngine.popen_uci("C:/path/to/stockfish.exe")
        # On Linux/macOS, it's often found if you install it via a package manager.
        engine = chess.engine.SimpleEngine.popen_uci("stockfish")
    except (FileNotFoundError, chess.engine.EngineTerminatedError):
        print("Stockfish engine not found.")
        print("Please download Stockfish from https://stockfishchess.org/download/")
        print("and ensure the executable is in your system's PATH or script directory.")
        return

    # Set up the board from the FEN string
    board = chess.Board(fen_string)

    # Analyze the position for 2 seconds to find the best move
    # The result object contains information about the engine's choice
    result = engine.play(board, chess.engine.Limit(time=2.0))
    best_move = result.move

    # Convert the move from UCI format (e.g., 'f2f7') to SAN format (e.g., 'Qxf7')
    san_move = board.san(best_move)

    print(f"For the position: {fen_string}")
    print(f"The best move for White is: {san_move}")

    # Close the engine
    engine.quit()


if __name__ == "__main__":
    # The FEN string for the given chess position
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    find_best_chess_move(fen)
