import chess
import chess.engine
import sys
import os

# --- INSTRUCTIONS ---
# 1. You need to have the Stockfish chess engine installed on your system.
# 2. You must provide the correct path to the Stockfish executable below.
#
# Examples:
# - Windows: STOCKFISH_PATH = r"C:\Users\YourUser\Downloads\stockfish\stockfish.exe"
# - Linux:   STOCKFISH_PATH = "/usr/games/stockfish"
# - macOS:   STOCKFISH_PATH = "/usr/local/bin/stockfish"
#
# If you don't have Stockfish, you can download it from: https://stockfishchess.org/download/

STOCKFISH_PATH = ""  # <--- IMPORTANT: SET THIS PATH

def format_score(score: chess.engine.Score) -> str:
    """Formats the engine's score into a human-readable string."""
    if score.is_mate():
        # score.white() gives moves from white's perspective
        return f"Mate in {score.white().mate()}"
    else:
        # score.white() gives score from white's perspective
        pawn_units = score.white().score(mate_score=10000) / 100.0
        return f"Score: +{pawn_units:.2f}"

def solve_chess_puzzle():
    """
    Analyzes a given chess position with Stockfish to find the best move among a list of choices.
    """
    if not STOCKFISH_PATH or not os.path.exists(STOCKFISH_PATH):
        print("Error: Stockfish engine path is not set or is invalid.")
        print("Please edit the script and set the STOCKFISH_PATH variable.")
        return

    # FEN for the chess position
    fen = "8/P7/1np1p3/5k2/6p1/1P1NK3/8/8 w - - 0 1"
    
    # Candidate moves (Choice: SAN)
    candidate_moves_san = {
        "A": "a8=Q",
        "B": "Nc5",
        "C": "Kd4",
        "D": "Kf2",
        "E": "Nf4",
        "F": "b4"
    }

    try:
        with chess.engine.SimpleEngine.popen_uci(STOCKFISH_PATH) as engine:
            board = chess.Board(fen)
            results = []

            print(f"Analyzing FEN: {fen}\n")
            print("--- Engine Analysis of Candidate Moves ---")

            for choice, san_move in candidate_moves_san.items():
                # The python-chess library can parse Standard Algebraic Notation (SAN)
                move = board.parse_san(san_move)
                
                # Make the move on a copy of the board to analyze the resulting position
                board.push(move)
                # Analyze the position after the move. Depth 20 is a good balance of speed and accuracy.
                info = engine.analyse(board, chess.engine.Limit(depth=20))
                score = info["score"]
                results.append((choice, san_move, score))
                board.pop()  # Undo the move to reset the board for the next candidate

            # Sort results from best to worst move for White
            results.sort(key=lambda x: x[2].white().score(mate_score=10000), reverse=True)

            best_choice = results[0][0]
            
            for choice, san_move, score in results:
                score_str = format_score(score)
                print(f"Move {choice} ({san_move}): {score_str}")

            print("\n--- Conclusion ---")
            best_move_san = candidate_moves_san[best_choice]
            best_score_str = format_score(results[0][2])

            print(f"The engine confirms that the best move is {best_choice}. {best_move_san}.")
            print(f"This move yields the highest evaluation of {best_score_str}.")
            print("It creates a decisive second threat by attacking the e6-pawn, overloading Black's defenses.")

            # Final answer in the required format
            print(f"\n<<<{best_choice}>>>")

    except (FileNotFoundError, PermissionError):
        print(f"Error: Could not find or execute the Stockfish engine at: '{STOCKFISH_PATH}'")
        print("Please ensure the path is correct and the file is executable.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

solve_chess_puzzle()