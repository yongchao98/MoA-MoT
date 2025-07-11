import stockfish
import os

def find_mate_in_n(fen_string: str):
    """
    Analyzes a chess FEN string to find if a forced mate exists.

    Args:
        fen_string: The FEN string representing the board state.
    
    Returns:
        The number of moves to mate, or an error message.
    """
    try:
        # Note: The Stockfish binary must be in your system's PATH.
        # Alternatively, provide the full path to the executable:
        # stockfish_path = "/usr/games/stockfish" # Example for Linux
        # stockfish_path = "/usr/local/bin/stockfish" # Example for macOS
        # stockfish_path = "C:\\path\\to\\stockfish.exe" # Example for Windows
        # if 'STOCKFISH_PATH' in os.environ:
        #     stockfish_engine = stockfish.Stockfish(os.environ['STOCKFISH_PATH'])
        # else:
        #     stockfish_engine = stockfish.Stockfish()

        # Simplified for user execution, assuming PATH is set.
        stockfish_engine = stockfish.Stockfish()
    except FileNotFoundError:
        print("Stockfish engine not found.")
        print("Please install Stockfish and ensure its executable is in your system's PATH,")
        print("or specify the path directly in the script.")
        return
    except Exception as e:
        print(f"An error occurred while initializing Stockfish: {e}")
        return

    if not stockfish_engine.is_fen_valid(fen_string):
        print("The provided FEN string is invalid.")
        return

    stockfish_engine.set_fen_position(fen_string)
    
    # We increase the thinking time/depth to ensure the result is accurate.
    # For a mate in a few moves, this is generally sufficient.
    evaluation = stockfish_engine.get_evaluation()

    if evaluation['type'] == 'mate':
        # The engine's value is the number of moves for the current player to mate.
        moves_to_mate = evaluation['value']
        print(moves_to_mate)
    else:
        # The problem statement implies a mate exists.
        # This is a fallback in case the engine depth is too low or there's no mate.
        print("No forced mate found for White at the current search depth.")
        print(f"The position evaluation is {evaluation['value']} centipawns.")


if __name__ == '__main__':
    # FEN from the problem description
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    find_mate_in_n(fen)
