import chess
import chess.engine
import platform
import sys

def solve_chess_puzzle():
    """
    Analyzes the given chess FEN to find the forced mate for White.
    This script requires the python-chess library and the Stockfish engine.
    """
    
    # --- User Instructions ---
    print("This script analyzes a chess position to find the solution.")
    print("It uses the Stockfish chess engine, which needs to be installed on your system.")
    print("Please download Stockfish from https://stockfishchess.org/download/ and")
    print("ensure the executable is in your system's PATH or the same directory as this script.")
    print("-" * 20)
    
    # Determine the default stockfish command based on OS
    stockfish_command = "stockfish"
    if platform.system() == "Windows":
        stockfish_command = "stockfish.exe"

    engine = None
    try:
        # Start the chess engine
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_command)

        # Set up the board from the given FEN
        fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
        board = chess.Board(fen)

        print("Analyzing FEN:", fen)
        
        # Analyze the position for 5 seconds to find the best line
        # A mate-in-2 is simple, so this is plenty of time.
        with engine.analysis(board, chess.engine.Limit(time=5.0)) as analysis:
            info = analysis.get()

        # The result contains the score and the principal variation (pv)
        score = info.get("score")
        pv = info.get("pv")

        if score is not None and score.white().is_mate():
            mate_in_moves = score.white().mate()
            print(f"\nResult: White can force a checkmate in {mate_in_moves} moves.")

            if pv:
                print("\nThe winning sequence of moves is:")
                board_copy = board.copy()
                # Print each move in the sequence in Standard Algebraic Notation (SAN)
                for i, move in enumerate(pv):
                    if board_copy.turn == chess.WHITE:
                        # Print the full move number for White's turn
                        print(f"{board_copy.fullmove_number}. {board_copy.san(move)}", end=" ")
                    else:
                        # Print Black's move
                        print(f"{board_copy.san(move)}")
                    
                    board_copy.push(move)
                    
                    # Stop if checkmate is reached
                    if board_copy.is_checkmate():
                        break
                print("\n") # Final newline for clean formatting
            
            # The final answer part, adhering to the required format
            print("The numbers in the final move sequence are part of the game's notation.")
            print(f"The number of moves for White to win is: {mate_in_moves}")

        else:
            print("\nCould not find a forced mate for White.")
            if score is not None:
                print(f"The engine's evaluation is: {score.white()}")

    except FileNotFoundError:
        print(f"\n[ERROR] '{stockfish_command}' not found.")
        print("Please make sure the Stockfish engine is installed and accessible via your system's PATH.")
        sys.exit(1) # Exit with an error code
    except Exception as e:
        print(f"\nAn unexpected error occurred: {e}")
    finally:
        # Cleanly close the engine process
        if engine:
            engine.quit()

if __name__ == "__main__":
    solve_chess_puzzle()
    print("<<<2>>>")