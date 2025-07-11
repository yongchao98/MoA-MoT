# This script requires the 'stockfish' python package and the Stockfish engine.
# 1. Install the package: pip install stockfish
# 2. Download Stockfish engine from https://stockfishchess.org/download/
#    Extract the archive and place the executable in a known location.

from stockfish import Stockfish
import os

def find_best_chess_move(fen_string):
    """
    Analyzes a chess position and finds the best move using the Stockfish engine.
    It will print an error message if the Stockfish executable is not found.
    """
    # --- Configuration ---
    # IMPORTANT: Update this path if Stockfish is not in your system's PATH
    # or the current directory.
    # Examples:
    # stockfish_path = "C:\\Users\\YourUser\\Downloads\\stockfish\\stockfish.exe"
    # stockfish_path = "/home/user/stockfish/stockfish"
    stockfish_path = "stockfish"

    # Add .exe for Windows if not present
    if os.name == 'nt' and not stockfish_path.endswith('.exe'):
        stockfish_path_exe = f"{stockfish_path}.exe"
        if os.path.exists(stockfish_path_exe):
            stockfish_path = stockfish_path_exe

    if not os.path.exists(stockfish_path):
        print(f"Error: Stockfish engine not found at '{stockfish_path}'.")
        print("Please download it from https://stockfishchess.org/download/")
        print("Then update the 'stockfish_path' variable in this script,")
        print("or place the executable in the same directory as the script.")
        return None

    try:
        # Initialize Stockfish with a reasonable depth for puzzle solving
        # Higher depth = stronger analysis, but more time.
        stockfish = Stockfish(path=stockfish_path, depth=18, parameters={"Threads": 2, "Hash": 1024})
        
        stockfish.set_fen_position(fen_string)
        
        # Get the best move from the engine
        best_move = stockfish.get_best_move()
        return best_move

    except Exception as e:
        print(f"An error occurred while using Stockfish: {e}")
        return None

# --- Main execution ---
if __name__ == "__main__":
    # The FEN string for the given chess position
    fen_position = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"

    print(f"Analyzing FEN: {fen_position}")
    
    # Find and print the best move
    best_move = find_best_chess_move(fen_position)
    
    if best_move:
        print(f"\nAssuming optimal play, White's best move is: {best_move}")
