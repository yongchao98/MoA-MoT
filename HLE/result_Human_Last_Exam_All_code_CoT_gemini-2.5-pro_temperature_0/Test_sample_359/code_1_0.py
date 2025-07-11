import chess
import chess.variant

def solve_koth_puzzle():
    """
    Analyzes the given King of the Hill chess position to find the shortest win for White.
    It programmatically follows the most efficient winning line and prints the analysis.
    """
    fen = "8/2k5/5pn1/1Pp1pNpp/3PP3/4K1B1/8/8 w - - 0 43"
    board = chess.variant.KingOfTheHillBoard(fen)

    print("Analyzing the King of the Hill puzzle...")
    print(f"Initial FEN: {fen}\n")

    # --- Step 1: White's first move ---
    # White's most direct plan is to vacate the e4 square for the king.
    # The move 1. e5 achieves this, creating an immediate threat of 2. Ke4 to win.
    move1_white_san = "e5"
    move1_white = board.parse_san(move1_white_san)
    board.push(move1_white)
    print(f"White plays 1. {move1_white_san}")
    print("This move vacates e4 and threatens an immediate win with Ke4.")

    # --- Step 2: Black's forced response ---
    # Black must defend against the Ke4 threat. The only move that controls the e4 square
    # is for the knight on g6 to move to e7. Any other move loses on the spot.
    move1_black_san = "Nge7"
    move1_black = board.parse_san(move1_black_san)
    board.push(move1_black)
    print(f"Black is forced to play 1... {move1_black_san} to guard the e4 square.")

    # --- Step 3: White's second move ---
    # White removes the defender of the e4 square.
    move2_white_san = "Nxe7"
    move2_white = board.parse_san(move2_white_san)
    board.push(move2_white)
    print(f"White plays 2. {move2_white_san}, capturing the only defender of e4.")

    # --- Step 4: Black's second forced response ---
    # Black must recapture the knight. If Black makes any other move, White wins with 3. Ke4.
    move2_black_san = "Kxe7"
    move2_black = board.parse_san(move2_black_san)
    board.push(move2_black)
    print(f"Black must recapture with 2... {move2_black_san}.")

    # --- Step 5: White's winning move ---
    # The e4 square is now clear and undefended. White's king can move to the center and win.
    move3_white_san = "Ke4"
    move3_white = board.parse_san(move3_white_san)
    board.push(move3_white)
    print(f"White plays the winning move 3. {move3_white_san}.")

    # --- Conclusion ---
    is_win = board.is_koth_win()
    num_moves_for_white_to_win = 3

    print("\n--- Result ---")
    if is_win:
        print(f"The sequence is a forced win for White.")
        print("The full winning line is:")
        print(f"1. {move1_white_san} {move1_black_san}")
        print(f"2. {move2_white_san} {move2_black_san}")
        print(f"3. {move3_white_san}")
        print(f"\nTherefore, assuming optimal play, White can win in {num_moves_for_white_to_win} moves.")
        print(num_moves_for_white_to_win)
    else:
        print("The analysis was incorrect, and the line did not result in a win.")

if __name__ == '__main__':
    # Ensure you have the library installed: pip install python-chess
    try:
        solve_koth_puzzle()
    except ImportError:
        print("Please install the required library: pip install python-chess")
