import chess
import chess.engine
import os

# --- Plan ---
# 1. Define the chess position using a FEN string.
# 2. IMPORTANT: Specify the correct path to your Stockfish engine executable below.
# 3. Use the engine to analyze the position and find the shortest mating sequence.
# 4. Print the mating sequence and the total number of White's moves to deliver mate.

# !!! IMPORTANT !!!
# Please replace the path below with the correct path to your Stockfish executable.
# For Linux/macOS, it might be "/usr/games/stockfish", "/usr/bin/stockfish" or "/usr/local/bin/stockfish".
# For Windows, it might be "C:/path/to/stockfish/stockfish.exe" (use forward slashes).
STOCKFISH_PATH = "/usr/games/stockfish"

# A helper to find stockfish if the default path doesn't work
def find_engine_or_fail(path):
    if os.path.exists(path):
        return path
    
    # Check common alternative paths if the primary one fails
    alternatives = [
        "/usr/bin/stockfish",
        "/usr/local/bin/stockfish",
        "stockfish" # In case it's in the system PATH
    ]
    for alt_path in alternatives:
        try:
            with chess.engine.SimpleEngine.popen_uci(alt_path) as engine:
                return alt_path
        except (FileNotFoundError, PermissionError):
            continue
    return None

def solve_checkmate():
    """
    Analyzes the chess position to find the number of moves to mate.
    """
    engine_path = find_engine_or_fail(STOCKFISH_PATH)
    if not engine_path:
        print(f"Error: Stockfish engine not found at '{STOCKFISH_PATH}' or in common locations.")
        print("Please install Stockfish and update the STOCKFISH_PATH variable in the script.")
        return

    # FEN string representing the position described in the problem.
    # r(a8)n(b8)__r(f8)k(g8) / p(a7)B(b7)p(c7)p(d7)Q(e7)_p(g7)p(h7) / 
    # _p(b6)__p(e6)b(f6)__ / ____N(e5)__Q(h5) / ___P(d4)N(e4)___ / __B(d3)_____ / 
    # PPP(a2,b2,c2)__PPP(f2,g2,h2) / R(a1)___K(e1)__R(h1) - White to move.
    board = chess.Board("rn2r1k1/pbppq1pp/1p2pb2/4N2Q/3PN3/2B5/PPP2PPP/R3K2R w KQ - 0 1")

    try:
        with chess.engine.SimpleEngine.popen_uci(engine_path) as engine:
            # Analyze the position to find the best line.
            # A depth of 20 is more than enough for a tactical puzzle.
            info = engine.analyse(board, chess.engine.Limit(depth=20))
            
            score = info.get("score").white()
            pv = info.get("pv") # Principal variation (the sequence of best moves)

            if score.is_mate():
                mate_in_moves = score.mate()
                
                print(f"White can deliver checkmate in {mate_in_moves} moves.")
                print("The final sequence of moves is:")
                
                # Print the moves in Standard Algebraic Notation (SAN)
                move_number = 1
                full_line_str = ""
                temp_board = board.copy()

                for move in pv:
                    if temp_board.turn == chess.WHITE:
                        # Add the move number for White's turn
                        full_line_str += f"{move_number}. {temp_board.san(move)} "
                    else:
                        full_line_str += f"{temp_board.san(move)} "
                        move_number += 1
                    
                    temp_board.push(move)
                    if temp_board.is_checkmate():
                        break # Stop printing moves after mate is reached
                
                # Add the checkmate symbol '#' to the end of the final move
                print(full_line_str.strip() + "#")

            else:
                print("No forced mate was found in the analysis.")

    except Exception as e:
        print(f"An error occurred while communicating with the chess engine: {e}")

# Run the solver
solve_checkmate()