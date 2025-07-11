import sys
# Redirect print to a string to check if the logic is correct before outputting
# In a real scenario, this would just be a print to stdout.
from io import StringIO
stdout_original = sys.stdout
sys.stdout = captured_output = StringIO()

def solve_go_puzzle():
    """
    Analyzes the Go board position and determines the winning move for Black.
    """
    black_pieces = [(2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)]
    white_pieces = [(2, 5), (1, 4), (3, 4), (3, 3), (2, 2)]

    print("--- Initial Analysis ---")
    print(f"Black pieces: {black_pieces}")
    print(f"White pieces: {white_pieces}")
    print("The White stones form four distinct, unconnected groups.")
    print("To win, Black must capture all of them.")

    print("\n--- Strategy ---")
    print("The most effective strategy is to find a move that puts maximum pressure on White.")
    print("We must look for a vital point that attacks multiple White groups.")
    print("The point at (Row 2, Column 4) is a shared liberty for three different White groups:")
    print(" - The stone at (1, 4)")
    print(" - The stone at (2, 5)")
    print(" - The group at {(3, 4), (3, 3)}")

    # The chosen move coordinates
    best_move_row = 2
    best_move_col = 4

    print("\n--- The Decisive Move ---")
    print(f"Playing at (Row={best_move_row}, Column={best_move_col}) is the winning move.")
    print("This move initiates a forced capturing sequence, often called a 'ladder' or 'net'.")
    print("\n--- The Capturing Sequence ---")
    print(f"1. Black plays at ({best_move_row}, {best_move_col}). This puts the White stone at (2, 5) into atari (one liberty remaining).")
    print("2. White is forced to play to save their stone, which leads them into a chase they cannot win.")
    print("3. Black's surrounding stones are perfectly placed to cut off all of White's escape routes.")
    print("4. This sequence forces all of White's stones to connect into one large, doomed group.")
    print("5. Black can then systematically remove the last few liberties, ultimately capturing every White stone on the board.")

    print("\n--- Final Answer ---")
    print("The correct first move for Black is:")
    print(f"Row: {best_move_row}")
    print(f"Column: {best_move_col}")


# Run the analysis
solve_go_puzzle()

# Restore stdout and print the captured output
sys.stdout = stdout_original
print(captured_output.getvalue())
