import sys

# This script solves the chess puzzle by outlining the forced checkmate sequence.
# It doesn't require a full chess engine, as the solution can be found through tactical analysis.

# Define a function to print the solution steps.
def solve_chess_puzzle():
    """
    Prints the step-by-step solution for the checkmate.
    """
    # The number of moves to mate
    num_moves = 2
    
    print(f"White can deliver checkmate in {num_moves} moves.")
    print("-" * 30)

    # Move 1
    print("Move 1 (White): Queen H5 x H7 check")
    print("This is a queen sacrifice that forces the Black king onto the h-file.")
    print("Move 1 (Black): King G8 x H7")
    print("Black's only legal move is to capture the queen.")
    print("-" * 30)

    # Move 2
    print("Move 2 (White): Knight E5 -> G6 checkmate")
    print("This move creates a double check from the knight and the rook on H1 (discovered check).")
    print("The Black king on H7 has no legal moves, resulting in checkmate.")
    print("-" * 30)

# Execute the function to display the solution.
if __name__ == "__main__":
    solve_chess_puzzle()
    # The final answer required by the user prompt
    sys.stdout.write("<<<2>>>\n")
