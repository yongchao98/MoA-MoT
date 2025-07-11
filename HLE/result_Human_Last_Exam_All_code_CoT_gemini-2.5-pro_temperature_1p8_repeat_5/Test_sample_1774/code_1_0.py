import chess
import chess.engine
import shutil
import sys

def find_stockfish_path():
    """Tries to find the stockfish executable in common locations."""
    # Check for stockfish in PATH
    path = shutil.which("stockfish")
    if path:
        return path

    # Check common paths for different OS
    possible_paths = [
        "./stockfish",  # In the current directory
        "/usr/games/stockfish",  # Common on Linux
        "/usr/bin/stockfish",
    ]
    for p in possible_paths:
        try:
            # Check if the file exists and is executable
            if chess.engine.os.access(p, chess.engine.os.X_OK):
                return p
        except Exception:
            continue
    return None

def solve_chess_puzzle():
    """
    Analyzes the given chess position to find the number of moves to mate.
    """
    # FEN representation of the board state based on the problem description.
    # Board: rn4kr/pbppq1pp/1p3b2/4N3/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1
    # White to move.
    fen_string = "rn4kr/pbppq1pp/1p3b2/4N3/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1"
    board = chess.Board(fen_string)

    stockfish_path = find_stockfish_path()
    if not stockfish_path:
        print("Stockfish engine not found.")
        print("Please install Stockfish (e.g., 'sudo apt install stockfish')")
        print("Or download it from https://stockfishchess.org/ and place the executable in this directory.")
        return

    try:
        with chess.engine.SimpleEngine.popen_uci(stockfish_path) as engine:
            # Analyze the position with a time limit to find the best move.
            # A 5-second analysis is more than sufficient for a tactical puzzle.
            info = engine.analyse(board, chess.engine.Limit(time=5.0))

            score = info.get("score")
            pv = info.get("pv")

            if score and pv and score.is_mate():
                # The score is from White's perspective.
                mate_in_moves = score.white().mate()
                print(f"White can force checkmate in {mate_in_moves} moves.")

                # Format and print the mating sequence.
                # The problem asked to "output each number in the final equation",
                # which is interpreted here as showing the full move notation.
                print("The winning sequence of moves is:")
                move_list = []
                # Make a copy of the board to play out the moves
                temp_board = board.copy()
                move_num = 1
                for i, move in enumerate(pv):
                    if temp_board.turn == chess.WHITE:
                        move_list.append(f"{move_num}. {temp_board.san(move)}")
                        move_num += 1
                    else:
                        # Append black's move to the last entry
                        move_list[-1] += f" {temp_board.san(move)}"
                    temp_board.push(move)
                
                print(" ".join(move_list))

            else:
                print("No forced checkmate sequence found in the given time.")
                if score:
                    print(f"Current evaluation: {score.white()}")

    except Exception as e:
        print(f"An error occurred while running the chess engine: {e}")
        print("Please ensure the Stockfish executable path is correct and it can be run.")

if __name__ == "__main__":
    solve_chess_puzzle()
