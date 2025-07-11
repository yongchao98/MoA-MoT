from stockfish import Stockfish

def solve_chess_mate():
    """
    Analyzes the given FEN string to find the number of moves for a forced win.

    This function requires the Stockfish chess engine to be installed and accessible
    on the system's PATH or at a specified location.
    """
    # You might need to change this path to where your Stockfish executable is located.
    # e.g., "C:\\Program Files\\stockfish\\stockfish.exe" on Windows or "/usr/bin/stockfish" on Linux.
    stockfish_path = "stockfish" 
    
    fen_position = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"

    try:
        # We increase the thinking time and depth to ensure the engine finds the mate.
        stockfish = Stockfish(path=stockfish_path, depth=20, parameters={"Threads": 2, "Minimum Thinking Time": 500})
        stockfish.set_fen_position(fen_position)
        
        evaluation = stockfish.get_evaluation()

        if evaluation['type'] == 'mate':
            moves_to_win = evaluation['value']
            # The prompt asks for each number in the final equation. 
            # The equation is simply that the number of moves to win is the found value.
            print(f"White can win in the following number of moves: {moves_to_win}")
        else:
            print("According to the engine, there is no forced mate from this position.")

    except (FileNotFoundError, PermissionError):
        print(f"Error: Stockfish executable not found or not accessible at the path: '{stockfish_path}'")
        print("Please install the Stockfish chess engine and ensure the path is correct.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

if __name__ == '__main__':
    solve_chess_mate()
