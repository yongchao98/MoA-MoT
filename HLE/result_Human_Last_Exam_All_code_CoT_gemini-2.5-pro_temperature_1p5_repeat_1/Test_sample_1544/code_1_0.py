# The user needs to have the stockfish engine installed and accessible in their system's PATH.
# For example, on Linux/macOS: sudo apt-get install stockfish or brew install stockfish
# On Windows, you may need to specify the path to the executable, e.g., Stockfish(path="/path/to/stockfish.exe")
try:
    from stockfish import Stockfish
except ImportError:
    print("The 'stockfish' library is not installed. Please install it using: pip install stockfish")
    exit()

def find_best_move():
    """
    Analyzes the given chess position to find the best move for White.
    """
    # The FEN from the prompt '8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1' is invalid.
    # The most plausible interpretation is a typo where '1kpP4' should be '1k1P4'.
    # We will use the corrected, valid FEN.
    fen_position = "8/3p4/1k1P4/p1q5/P7/8/5Q2/6K1 w - - 0 1"

    try:
        # Initialize Stockfish. Adjust the path if the executable is not in your system's PATH.
        stockfish = Stockfish(depth=18, parameters={"Threads": 2, "Minimum Thinking Time": 1000})
        
        stockfish.set_fen_position(fen_position)
        
        best_move = stockfish.get_best_move()
        
        print(f"The FEN of the position is: {fen_position}")
        print("Assuming optimal play from both sides...")
        print(f"White's best move is: {best_move}")

    except Exception as e:
        print(f"An error occurred. Please ensure the Stockfish engine is installed and accessible.")
        print(f"Error details: {e}")

if __name__ == "__main__":
    find_best_move()