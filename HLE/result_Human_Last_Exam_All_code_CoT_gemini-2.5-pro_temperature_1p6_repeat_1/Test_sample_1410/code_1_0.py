import chess
import chess.engine
import sys
import os

def solve_chess_puzzle(fen_string, engine_path):
    """
    Analyzes a chess position from a FEN string to find the mate-in-N solution.

    Args:
        fen_string (str): The FEN string representing the board state.
        engine_path (str): The file path to the UCI chess engine executable.
    """
    if not os.path.exists(engine_path):
        print(f"Error: Chess engine not found at '{engine_path}'")
        print("Please download a UCI engine (like Stockfish) and update the 'engine_path'.")
        return

    try:
        engine = chess.engine.SimpleEngine.popen_uci(engine_path)
        board = chess.Board(fen_string)

        # Analyze the position to find the mate
        # A time limit of 5 seconds is more than enough for most mate puzzles.
        info = engine.analyse(board, chess.engine.Limit(time=5.0))

        score = info.get("score")
        if score is None or not score.is_mate():
            print("Could not find a forced mate in the given time.")
            return

        # The mate score is given from the perspective of the side to move.
        mate_in = score.white().mate()
        
        # Get the principal variation (the sequence of best moves)
        pv = info.get("pv")
        if pv:
            print("The winning move sequence is:")
            # Convert moves to standard algebraic notation for readability
            move_sequence = " ".join([board.san(move) for move in pv])
            # The mate is delivered on the Nth move by white. Let's show the full sequence leading to it.
            moves_to_show = mate_in * 2 -1
            print(f"1. {move_sequence.split(' ')[0]} {move_sequence.split(' ')[1]} 2. {move_sequence.split(' ')[2]} {move_sequence.split(' ')[3]} 3. {move_sequence.split(' ')[4]}")


        print("\n-------------------------------------------")
        print(f"White can force a checkmate in {mate_in} moves.")
        print("-------------------------------------------")


    except Exception as e:
        print(f"An error occurred: {e}")
    finally:
        if 'engine' in locals():
            engine.quit()

if __name__ == "__main__":
    # The FEN string from the problem
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    
    # --- IMPORTANT ---
    # You must specify the path to your Stockfish (or other UCI engine) executable here.
    # Examples: "/usr/games/stockfish", "C:/engines/stockfish.exe", "./stockfish"
    # On linux you can find it by typing 'which stockfish' in your terminal
    if sys.platform == "win32":
        # Example for Windows
        engine_path = "stockfish.exe"
    else:
        # Example for Linux/macOS
        engine_path = "/usr/games/stockfish" 

    # Attempt to find stockfish automatically if the user-provided path doesn't exist
    if not os.path.exists(engine_path):
        try:
            import shutil
            found_path = shutil.which("stockfish")
            if found_path:
                engine_path = found_path
            else: # last attempt for common location
                 if os.path.exists("/usr/bin/stockfish"):
                     engine_path = "/usr/bin/stockfish"
        except ImportError:
            pass
            
    solve_chess_puzzle(fen, engine_path)
