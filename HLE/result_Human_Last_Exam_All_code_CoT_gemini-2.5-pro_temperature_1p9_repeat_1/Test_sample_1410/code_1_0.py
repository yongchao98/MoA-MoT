import os
import sys
from stockfish import Stockfish

def find_mate_in_n():
    """
    Analyzes a chess FEN to find if a forced checkmate exists
    and in how many moves.
    """
    # This script requires the Stockfish engine.
    # We will try a few common executable names.
    # Please install Stockfish and, if necessary, update the path below.
    stockfish_path = None
    possible_paths = ["stockfish", "/usr/bin/stockfish", "/usr/local/bin/stockfish"]
    if sys.platform == "win32":
        # Add common windows paths/names if necessary
        possible_paths.extend(["stockfish.exe", "stockfish_20090216_x64.exe"])

    for path in possible_paths:
        # Using shutil.which to find the executable in the system's PATH
        try:
            import shutil
            if shutil.which(path):
                stockfish_path = path
                break
        except ImportError:
            # Fallback for systems without shutil
            if os.path.exists(path) and os.access(path, os.X_OK):
                 stockfish_path = path
                 break
    
    if not stockfish_path:
        print("Error: Stockfish executable not found.")
        print("Please install Stockfish from https://stockfishchess.org/download/")
        print("and ensure its location is in your system's PATH, or update the 'possible_paths' in this script.")
        return

    try:
        # We set a sufficient depth for finding mates in typical chess puzzles.
        stockfish = Stockfish(path=stockfish_path, depth=20, parameters={"Threads": 2, "Hash": 128})
    except (FileNotFoundError, PermissionError) as e:
        print(f"Error initializing Stockfish: {e}")
        return

    fen_string = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    if not stockfish.is_fen_valid(fen_string):
        print(f"The provided FEN string is invalid: {fen_string}")
        return
        
    stockfish.set_fen_position(fen_string)

    # Get the evaluation from the engine
    evaluation = stockfish.get_evaluation()

    print(f"Analyzing chess position: {fen_string}")
    print("Assuming optimal play from both White and Black...")

    if evaluation['type'] == 'mate':
        moves_to_mate = evaluation['value']
        if moves_to_mate > 0:
            print(f"\nWhite has a forced checkmate in {moves_to_mate} moves.")
            # The puzzle asks for the number of moves, which is our result.
            # There is no "equation", so we just output the number itself.
            print("\nFinal Answer (Number of moves for White to win):")
            print(moves_to_mate)

        else:
            # This case means Black is the one with the mate
            print(f"\nBlack has a forced checkmate in {-moves_to_mate} moves.")
    else:
        # In case no mate is found, we show the pawn advantage evaluation
        cp_advantage = evaluation['value'] / 100.0
        print(f"\nNo forced mate was found. The position evaluation is {cp_advantage} in favor of White.")

# Run the analysis
find_mate_in_n()