import chess
import chess.engine
import sys
import os

def solve_chess_puzzle():
    """
    Sets up a chess position and uses an engine to find the mate in X moves.
    """
    # Define the path to the stockfish executable.
    # IMPORTANT: Change this to the actual path of your stockfish executable,
    # or ensure 'stockfish' is in your system's PATH.
    # For example on Windows: "C:\\Users\\YourUser\\Downloads\\stockfish\\stockfish.exe"
    # For example on Linux/Mac: "/home/user/stockfish/stockfish"
    stockfish_executable = "stockfish" 

    # Verify if the stockfish executable exists at the given path or in the system PATH
    def is_in_path(program):
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                return True
        return False

    if not (os.path.isfile(stockfish_executable) or is_in_path(stockfish_executable.split(os.sep)[-1])):
        print(f"Error: Stockfish engine not found at '{stockfish_executable}'.", file=sys.stderr)
        print("Please download Stockfish from https://stockfishchess.org/download/", file=sys.stderr)
        print("Then, update the 'stockfish_executable' variable in the script with the correct file path.", file=sys.stderr)
        return

    # FEN for the position: rn3rk1/pbpq1pp1/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1
    # This string represents the exact board state you described.
    board = chess.Board("rn3rk1/pbpq1pp1/1p2pb2/4N2Q/3PN3/3B4/PPP2PPP/R3K2R w KQ - 0 1")

    engine = None
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_executable)

        # Analyze the position to find the best move and check for a mate sequence.
        # We'll give the engine a few seconds to ensure it finds the mate.
        info = engine.analyse(board, chess.engine.Limit(time=3.0))

        # The engine provides a score. If it's a mate, the score will be a Mate object.
        score = info.get("score").white()
        
        if not score.is_mate():
            print("No forced checkmate found in a reasonable time.")
            return

        mate_in = score.mate()
        # The engine gives the number of moves until mate, which includes black's responses.
        # So a "mate in 2" means W-B-W(mate). We report this as 2 moves.
        white_moves_to_mate = mate_in

        print(f"White has a forced checkmate in {white_moves_to_mate} moves.")
        print("The winning equation is:")

        # The 'pv' (principal variation) is the sequence of best moves.
        best_move_sequence = info.get("pv")

        move_counter = 1
        # We iterate through the sequence of moves and print them in order.
        for i in range(len(best_move_sequence)):
            move = best_move_sequence[i]
            # Print White's move with the move number
            if board.turn == chess.WHITE:
                print(f"{move_counter}. {board.san(move)}", end=" ")
                move_counter += 1
            # Print Black's response
            else:
                print(f"{board.san(move)}", end="\n")
            
            board.push(move)

        # Ensure the last line ends with a newline
        print("\n")
        print("Final Position:")
        print(board)
        
        # The final answer as requested
        print(f"<<<{white_moves_to_mate}>>>")

    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
    finally:
        if engine:
            engine.quit()

if __name__ == "__main__":
    solve_chess_puzzle()