import chess
import chess.engine
import sys

# --- Configuration ---
# IMPORTANT: Please update this path to the location of your Stockfish engine executable.
# You can download Stockfish from https://stockfishchess.org/download/
# For Linux/macOS, it might be '/usr/games/stockfish' or a local path like './stockfish'.
# For Windows, it might be 'C:/path/to/stockfish.exe'.
STOCKFISH_PATH = "/usr/games/stockfish" # <-- CHANGE THIS IF NEEDED

def check_stockfish_path():
    """Checks if the Stockfish executable exists at the specified path."""
    import os
    if not os.path.exists(STOCKFISH_PATH):
        print(f"Error: Stockfish engine not found at '{STOCKFISH_PATH}'.")
        print("Please download Stockfish from https://stockfishchess.org/download/")
        print("and update the STOCKFISH_PATH variable in the script.")
        sys.exit(1)

def evaluate_best_move():
    """
    Analyzes a given chess position with Stockfish to find the best move.
    """
    check_stockfish_path()
    
    try:
        engine = chess.engine.SimpleEngine.popen_uci(STOCKFISH_PATH)
    except FileNotFoundError:
        print(f"Error: Could not start the engine at '{STOCKFISH_PATH}'.")
        print("Please ensure the path is correct and the file is executable.")
        return
        
    # Setup the board from the FEN string
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    board = chess.Board(fen)

    moves_san = {
        "A": "a8=Q",
        "B": "Nc5",
        "C": "Kd4",
        "D": "Kf2",
        "E": "Nf4",
        "F": "b4",
    }
    
    print("Analyzing White's best move...")
    print("Position FEN:", fen)
    print("-" * 30)

    best_move_label = None
    best_score = None

    for label, san_move in moves_san.items():
        try:
            # Convert move from Standard Algebraic Notation (SAN) to a move object
            move = board.parse_san(san_move)
            
            # Analyze the position after this move is made
            info = engine.analyse(board, chess.engine.Limit(time=0.5), root_moves=[move])
            score = info["score"].white()
            
            # Print the evaluation for the move
            # The score is from White's perspective.
            # PovScore(Cp(+650), WHITE) means White is up by 6.5 pawns.
            # PovScore(Mate(+5), WHITE) means White has a mate in 5.
            print(f"Move {label} ({san_move}): Evaluation = {score}")

            # Keep track of the best move found so far
            if best_score is None or score > best_score:
                best_score = score
                best_move_label = label

        except chess.InvalidMoveError:
            print(f"Move {label} ({san_move}): Is an invalid move in this position.")
        except Exception as e:
            print(f"An error occurred while analyzing move {san_move}: {e}")


    engine.quit()
    
    print("-" * 30)
    if best_move_label:
        print(f"The best move among the choices is {best_move_label} ({moves_san[best_move_label]}) with a score of {best_score}.")
    else:
        print("Could not determine the best move.")

    # Final answer in the required format
    # print(f"<<<{best_move_label}>>>")


if __name__ == "__main__":
    evaluate_best_move()
<<<B>>>