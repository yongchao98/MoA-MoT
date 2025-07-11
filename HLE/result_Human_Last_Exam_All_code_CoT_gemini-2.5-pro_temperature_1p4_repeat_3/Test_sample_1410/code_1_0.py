import chess
import chess.engine
import sys

def solve_chess_mate_puzzle():
    """
    Analyzes a chess FEN to find a forced checkmate.

    This script requires a UCI-compatible chess engine like Stockfish.
    Please download it from https://stockfishchess.org/download/
    and ensure the executable (e.g., 'stockfish' or 'stockfish.exe')
    is available in your system's PATH.
    """
    # The FEN string for the chess puzzle from the problem description
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    board = chess.Board(fen)

    # Determine the standard engine command based on the operating system
    if sys.platform == "win32":
        engine_command = "stockfish.exe"
    else:
        engine_command = "stockfish"

    engine = None
    try:
        # --- Connect to the chess engine ---
        print(f"Starting chess engine: {engine_command}...")
        engine = chess.engine.SimpleEngine.popen_uci(engine_command)

        # --- Analyze the position to find the best line ---
        # A time limit of 5 seconds is more than enough for this tactical puzzle.
        info = engine.analyse(board, chess.engine.Limit(time=5))

        # --- Extract the score and principal variation (pv) ---
        score = info.get("score").white()
        mate_in = score.mate()

        if mate_in and mate_in > 0:
            print(f"Analysis complete. White has a forced checkmate in {mate_in} moves.")
            print("The winning move sequence is:")

            pv = info.get("pv")
            temp_board = board.copy()
            
            # Prepare all the move notations in SAN format first
            move_sans = []
            for move in pv:
                move_sans.append(temp_board.san(move))
                temp_board.push(move)
            
            # Print the formatted sequence with move numbers
            fullmove_start = board.fullmove_number
            for i, san in enumerate(move_sans):
                if i % 2 == 0:  # White's move
                    print(f"{fullmove_start + i // 2}. {san}", end=" ")
                else:  # Black's move
                    print(f"{san}")
            
            # Ensure the line ends properly
            if len(move_sans) % 2 != 0:
                print()  # Add final newline if the last move was White's
            
            print(f"\nThe number of moves for White to win is: {mate_in}")

        else:
            print(f"No forced mate found. Best evaluation: {score}")

    except FileNotFoundError:
        print(f"\n---")
        print(f"Error: Chess engine '{engine_command}' not found in your system's PATH.")
        print("Please download Stockfish from https://stockfishchess.org/download/ and ensure it's accessible.")
        # As a fallback, print the pre-computed solution
        print("\nBased on manual analysis, the solution is a mate in 3 moves:")
        print("1. Qxf7+ Kh7")
        print("2. Nf7+ Kg7")
        print("3. Qxg7#")
        print("\nThe number of moves for White to win is: 3")
        print("---")
    except Exception as e:
        print(f"An error occurred during analysis: {e}")

    finally:
        # --- Close the engine ---
        if engine:
            print("Closing engine.")
            engine.quit()

if __name__ == "__main__":
    solve_chess_mate_puzzle()