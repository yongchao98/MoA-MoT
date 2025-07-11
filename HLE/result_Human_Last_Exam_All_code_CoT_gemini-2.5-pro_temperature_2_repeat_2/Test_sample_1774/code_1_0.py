import chess
import chess.engine
import os

def solve_chess_puzzle():
    """
    Analyzes the specified chess position to find the forced mate sequence.

    This script requires a UCI-compatible chess engine like Stockfish.
    Please download it from https://stockfishchess.org/download/ and
    place the executable file (e.g., "stockfish.exe") in the same folder
    as this Python script, or provide the full path to it below.
    """
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    # CONFIGURE ENGINE PATH:
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---
    # Try to find the engine executable automatically.
    # If the script fails to find it, update the `engine_path`
    # variable with the correct path to your Stockfish executable.
    # Example for Windows: r"C:\Users\YourUser\Downloads\stockfish\stockfish.exe"
    # Example for Linux/Mac: "/usr/games/stockfish" or "./stockfish"
    
    engine_path = "stockfish" # Assumes 'stockfish' is in your system's PATH
    
    # Simple check for the executable in the current directory if not in PATH
    if not os.path.exists(engine_path) or not os.access(engine_path, os.X_OK):
      local_paths = ["./stockfish", "./stockfish.exe", "stockfish"]
      found = False
      for path in local_paths:
          if os.path.exists(path) and os.access(path, os.X_OK):
              engine_path = path
              found = True
              break
      if not found:
          print("--- chess engine not found ---")
          print("Error: Could not find the Stockfish engine.")
          print("Please download Stockfish and place its executable in this directory or update the 'engine_path' variable in the script.")
          return
    # --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- --- ---

    try:
        engine = chess.engine.SimpleEngine.uci(engine_path)
    except FileNotFoundError:
        print("Error: The chess engine was not found at the specified path.")
        return
    except Exception as e:
        print(f"An error occurred while starting the engine: {e}")
        return

    # FEN string representing the board state from the problem description
    fen = "rn2r1k/pbppq1pp/1p2pb2/4N2Q/3PN3/3B4/PPP1P1PP/R3K2R w KQq - 0 1"
    board = chess.Board(fen)

    try:
        # Analyze the position to find the best line (limit search to 5 seconds)
        info = engine.analyse(board, chess.engine.Limit(time=5))
    except Exception as e:
        print(f"An error occurred during analysis: {e}")
        engine.quit()
        return
        
    pv = info.get("pv")
    score = info.get("score")

    if not pv or not score:
        print("Engine did not return a solution. Try increasing the search time.")
        engine.quit()
        return
        
    # Get the score from White's perspective
    white_score = score.white()

    if white_score.is_mate():
        mate_in = white_score.mate()
        
        # Build the final equation string
        final_equation = []
        # Copy the board to play out the moves without changing the original
        temp_board = board.copy()
        
        for i, move in enumerate(pv):
            # Format White's move with the move number
            if temp_board.turn == chess.WHITE:
                final_equation.append(f"{temp_board.fullmove_number}. {temp_board.san(move)}")
            # Format Black's move
            else:
                final_equation.append(temp_board.san(move))
            
            temp_board.push(move)

        # Add the checkmate symbol '#' to the very last move
        final_equation[-1] += "#"
        
        print(f"White has a forced checkmate in {mate_in} moves.")
        print("The final mating equation is:")
        print(" ".join(final_equation))
    else:
        best_move = board.san(pv[0])
        print(f"No forced mate was found. The best move is {best_move} with an evaluation of {white_score}.")

    engine.quit()

if __name__ == '__main__':
    solve_chess_puzzle()