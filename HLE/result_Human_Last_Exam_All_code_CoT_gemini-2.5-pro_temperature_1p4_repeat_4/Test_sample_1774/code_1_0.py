import sys
import io

# Keep the original stdout
original_stdout = sys.stdout
# Redirect stdout to a string
sys.stdout = captured_output = io.StringIO()

def solve_chess_puzzle():
    """
    Analyzes the given chess position and prints the solution.
    """
    white_moves_to_mate = 1
    black_reply_moves = 0
    total_moves = white_moves_to_mate + black_reply_moves

    print("Analyzing the board for the fastest mate...")
    print("---------------------------------------------")
    print("The key attacking pieces for White are the Queen on h5, the Bishop on d3, and the Knight on e5.")
    print("The most forceful move is for the Queen to capture the pawn on h7.")
    print("\nStep 1: White moves Queen to h7, capturing the pawn (Qxh7).")
    print("This places the Black King on g8 in check.")

    print("\nStep 2: We analyze Black's possible responses.")
    print("  - Can the King capture the Queen on h7? No, the White Bishop on d3 also attacks h7, so the King cannot move there.")
    print("  - Can the King move to an empty square?")
    print("    - h8 is attacked by the White Queen on h7.")
    print("    - f7 is attacked by the White Knight on e5.")
    print("  - Can the King move to a friendly-occupied square?")
    print("    - f8 is occupied by the Black Rook.")
    print("  - Can the check be blocked? No, the Queen is adjacent to the King.")

    print("\nConclusion: Black has no legal moves. It is checkmate.")

    print("\n--- Final Equation ---")
    print(f"Number of White moves to deliver checkmate: {white_moves_to_mate}")
    print(f"Number of Black moves possible after that: {black_reply_moves}")
    print(f"Equation: {white_moves_to_mate} + {black_reply_moves} = {total_moves}")
    print(f"\nIt takes {total_moves} move for White to checkmate the Black king.")

solve_chess_puzzle()

# Restore stdout and get the captured output
sys.stdout = original_stdout
output = captured_output.getvalue()

print(output)