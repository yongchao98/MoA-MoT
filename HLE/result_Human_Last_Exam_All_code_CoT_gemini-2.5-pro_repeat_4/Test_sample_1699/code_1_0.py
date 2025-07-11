import sys
import io

# Redirect stdout to capture the output
old_stdout = sys.stdout
sys.stdout = captured_output = io.StringIO()

def solve_go_problem():
    """
    This function analyzes the Go board position and determines the optimal move for Black.
    """
    # Step 1 & 2: Define the board state and analyze White's groups.
    # Black pieces: (2, 6), (4, 6), (3, 5), (4, 4), (4, 3), (4, 2), (6, 2)
    # White pieces: (2, 5), (1, 4), (3, 4), (3, 3), (2, 2)
    #
    # The White stones are separated into four small, vulnerable groups:
    # W1: {(2, 5)}
    # W2: {(1, 4)}
    # W3: {(3, 4), (3, 3)}
    # W4: {(2, 2)}
    #
    # To win, Black must capture all four groups. This requires preventing them
    # from connecting and forming a stable, "living" shape.

    # Step 3: Analyze Candidate Moves.
    # The most critical points on the board are those that affect multiple groups.
    # The point (2, 4) is adjacent to stones from groups W1, W2, and W3, making it a vital strategic point.
    # Let's analyze the sequence of plays starting with Black playing at (2, 4).

    print("Analyzing the consequences of Black playing at (2, 4):")
    print("-" * 50)

    # Step 4: Simulate the game after the key move.
    print("1. Black plays at (2, 4).")
    print("   - This move puts the White group W1 {(2, 5)} into 'atari' (one liberty remaining at (1, 5)).")
    print("   - White is forced to respond to save this stone.")
    print("\n2. White plays at (1, 5).")
    print("   - This saves W1 and connects it with W2, forming a new group G_A = {(1, 4), (1, 5), (2, 5)}.")
    print("   - G_A has two liberties: (1, 3) and (1, 6).")
    print("\n3. Black plays at (3, 2).")
    print("   - This is another tactical move, putting the White group W3 {(3, 4), (3, 3)} into atari.")
    print("   - White is again forced to respond.")
    print("\n4. White plays at (2, 3).")
    print("   - This saves W3 and connects it to W4, forming group G_B = {(2, 2), (2, 3), (3, 3), (3, 4)}.")
    print("   - G_B now has three liberties: (1, 2), (1, 3), and (2, 1).")
    print("\n5. Black plays the decisive move at (1, 3).")
    print("   - This single move attacks both of White's newly formed groups, as (1, 3) is a liberty for both G_A and G_B.")
    print("   - G_A is now in atari with only one liberty at (1, 6).")
    print("   - G_B is also in atari with liberties at (1, 2) and (2, 1).")
    print("\n6. White is in a 'miai' situation and cannot save both groups.")
    print("   - If White saves G_A by playing at (1, 6), Black will capture G_B.")
    print("   - If White tries to save G_B, Black will capture G_A.")
    print("   - After one group is captured, the other is easily captured as well.")

    # Step 5: Conclusion
    # The initial move at (2, 4) starts a chain reaction that inevitably leads to the capture of all white stones.
    final_move_row = 2
    final_move_col = 4

    print("-" * 50)
    print("Conclusion: The first move that guarantees victory is (2, 4).")
    # As per the instruction to "output each number in the final equation"
    print(f"Final Move Row = {final_move_row}")
    print(f"Final Move Column = {final_move_col}")

solve_go_problem()

# Restore stdout and print the captured output
sys.stdout = old_stdout
print(captured_output.getvalue())