import chess
import chess.engine
import sys
import os

def find_shortest_mate():
    """
    Analyzes a chess position from a FEN string to find the shortest
    forced checkmate.

    This script requires the Stockfish chess engine to be installed.
    It will try to find Stockfish in common locations, but you may
    need to set the STOCKFISH_PATH environment variable if it's not found.
    """
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)

    # --- Find Stockfish Engine ---
    # You can set the STOCKFISH_PATH environment variable to the path of the
    # Stockfish executable. Otherwise, the script will check common paths.
    stockfish_path = os.getenv("STOCKFISH_PATH")
    if not stockfish_path:
        # Common paths for different OS
        paths = [
            "/usr/local/bin/stockfish",
            "/usr/bin/stockfish",
            "/opt/homebrew/bin/stockfish", # For macOS with Homebrew
            "stockfish.exe", # For Windows if in the same folder or in PATH
            "stockfish" # For Linux if in PATH
        ]
        for path in paths:
            if os.path.exists(path) and os.access(path, os.X_OK):
                stockfish_path = path
                break

    if not stockfish_path:
        print("---")
        print("Error: Stockfish engine not found.")
        print("Please install Stockfish from https://stockfishchess.org/download/")
        print("and ensure it is in your system's PATH, or set the STOCKFISH_PATH environment variable.")
        print("---")
        sys.exit(1)
        
    engine = None
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)

        # Analyze the position. A short time limit is sufficient for tactical puzzles.
        analysis = engine.analyse(board, chess.engine.Limit(time=2.0))
        
        score = analysis.get("score").white()

        if score.is_mate():
            moves_to_mate = score.mate()
            winning_variation = analysis.get("pv")
            
            # To be clear, "mate in N" means N moves for the winning side.
            # We will display the full move sequence for clarity.
            print("Position analysis complete.")
            print(f"White has a forced checkmate in {moves_to_mate} moves.")
            
            print("\nThe winning sequence (Principal Variation) is:")
            # Display the move sequence in Standard Algebraic Notation (SAN)
            temp_board = board.copy()
            move_san_list = []
            for i, move in enumerate(winning_variation):
                if temp_board.turn == chess.WHITE:
                    move_number = temp_board.fullmove_number
                    move_san_list.append(f"{move_number}. {temp_board.san(move)}")
                else:
                    move_san_list[-1] += f" {temp_board.san(move)}"
                temp_board.push(move)
                # Stop printing after the checkmate
                if temp_board.is_checkmate():
                    break

            print(" ".join(move_san_list))
            
            # Fulfilling the "output each number in the final equation" instruction
            print("\nFinal Answer Calculation:")
            print(f"Number of moves for White to win = {moves_to_mate}")

        else:
            print(f"No forced mate was found. The position evaluation is: {score}")

    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if engine:
            engine.quit()

if __name__ == "__main__":
    find_shortest_mate()