import os
import shutil
from stockfish import Stockfish

def solve_chess_mate_puzzle():
    """
    This script solves for the number of moves to a forced win in a given chess position.
    
    It requires:
    1. The 'stockfish' python library: pip install stockfish
    2. The Stockfish chess engine executable installed and in the system's PATH,
       or in a common location like /usr/bin/stockfish.
    """
    
    fen_string = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    
    #--- Find Stockfish executable ---
    stockfish_path = None
    # Check common paths for Linux/macOS
    common_paths = [
        "/usr/games/stockfish",
        "/usr/bin/stockfish",
        "/opt/homebrew/bin/stockfish",
        "/usr/local/bin/stockfish"
    ]
    for path in common_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            stockfish_path = path
            break
            
    # Fallback to checking the system PATH (works for most installations, including Windows)
    if not stockfish_path:
        if shutil.which("stockfish"):
            stockfish_path = "stockfish"

    if not stockfish_path:
        print("Error: Could not find the Stockfish executable.")
        print("Please install Stockfish from https://stockfishchess.org/download/")
        print("and ensure it is in your system's PATH.")
        return

    try:
        # Initialize the engine with a sufficient search depth to find the mate.
        stockfish = Stockfish(path=stockfish_path, depth=22, parameters={"Threads": 2, "Minimum Thinking Time": 500})
        
        stockfish.set_fen_position(fen_string)
        
        evaluation = stockfish.get_evaluation()
        
        if evaluation.get('type') == 'mate':
            # A positive value means the current player (White) can force a mate.
            mate_in_moves = evaluation.get('value')
            if mate_in_moves > 0:
                # The question asks for the number of moves for White to win. This is the value.
                # The instruction about an "equation" is interpreted as printing the final number.
                print(mate_in_moves)
            else:
                print("The engine found a mate, but not for White.")
        else:
            print("No forced mate was found by the engine at this depth.")
            print("Position evaluation:", evaluation)
            
    except Exception as e:
        print(f"An error occurred while running Stockfish: {e}")

if __name__ == '__main__':
    solve_chess_mate_puzzle()
