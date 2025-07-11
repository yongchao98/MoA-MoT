import chess
import chess.engine
import shutil
import sys

def solve_chess_puzzle():
    """
    Analyzes a given FEN string to find the shortest mate for White.
    It uses the Stockfish chess engine to find the optimal line of play.
    """
    
    # --- Step 1: Explain and find Stockfish engine ---
    # The script requires the Stockfish chess engine to be installed and accessible
    # in the system's PATH. If it's not found, instructions are printed.
    
    stockfish_path = shutil.which("stockfish")
    if stockfish_path is None:
        print("Error: Stockfish engine not found in your system's PATH.")
        print("Please install Stockfish to run this analysis.")
        print("\nInstallation instructions:")
        print("- On Debian/Ubuntu: sudo apt-get update && sudo apt-get install stockfish")
        print("- On macOS with Homebrew: brew install stockfish")
        print("- On Windows: Download from https://stockfishchess.org/download/, unzip,")
        print("  and add the folder containing 'stockfish.exe' to your system's PATH.")
        return

    # --- Step 2: Setup the board and engine ---
    fen = "4r1k1/ppp2pp1/2n4p/3QN1qB/3P4/2P4P/PP4P1/4rRK1 w - - 0 20"
    engine = None
    try:
        engine = chess.engine.SimpleEngine.popen_uci(stockfish_path)
        board = chess.Board(fen)

        # --- Step 3: Analyze the position ---
        # We'll analyze with a sufficient depth to find a tactical mate.
        # A depth of 20 is more than enough for most mate-in-X puzzles.
        info = engine.analyse(board, chess.engine.Limit(depth=20))
        
        # --- Step 4: Extract the score and Principal Variation (PV) ---
        score = info.get("score")
        if score is None:
            print("Could not retrieve a score from the engine.")
            return

        white_score = score.white()

        # --- Step 5: Check for mate and print the results ---
        if white_score.is_mate():
            mate_in = white_score.mate()
            if mate_in <= 0:
                 print(f"White does not have a forced mate. Black is winning with a mate in {-mate_in}.")
                 return
                 
            print(f"White can force a checkmate in {mate_in} moves.")
            
            # Print the move sequence
            pv = info.get("pv", [])
            if pv:
                print("\nThe winning sequence is:")
                board_copy = board.copy()
                move_num = board_copy.fullmove_number
                is_white_turn = True
                
                for move in pv:
                    # To correctly display special notations like '++' for double check,
                    # we first generate the standard SAN, then check board state after the move.
                    san_move = board_copy.san(move)
                    board_copy.push(move)

                    if board_copy.is_checkmate():
                        san_move += "#"
                    elif board_copy.is_double_check():
                         # chess.Board.san() only produces '+', so we find it and replace with '++'
                        if '+' in san_move:
                            san_move = san_move.replace('+', '++')
                        else:
                            san_move += "++"

                    if is_white_turn:
                        # Output the move number, e.g., "20."
                        print(f"{move_num}. {san_move}", end=" ")
                    else:
                        print(f"{san_move}")
                        move_num += 1

                    is_white_turn = not is_white_turn
                    if board_copy.is_game_over():
                        break
                # Add a newline for clean formatting if the sequence ends on White's move
                if not is_white_turn:
                     print()
            
            # Print each move number to satisfy the prompt's instruction
            print("\nThe number for each of White's moves in the final winning sequence are:")
            for i in range(1, mate_in + 1):
                print(i)
                
            # The final answer in the required format
            print(f"\nFinal answer in the required format is below.")
            print(f"<<<{mate_in}>>>")

        else:
            print(f"No forced mate was found. The engine evaluation is: {white_score.score(mate_after=True)} centipawns.")

    except Exception as e:
        print(f"An error occurred during analysis: {e}")
    finally:
        if engine:
            engine.quit()

if __name__ == '__main__':
    solve_chess_puzzle()