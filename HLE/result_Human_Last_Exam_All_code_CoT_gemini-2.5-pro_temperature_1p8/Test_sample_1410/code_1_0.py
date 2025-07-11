import sys

# This script prints the solution to the chess puzzle.
# The puzzle is to find the minimum number of moves for White to win,
# which means finding the shortest forced checkmate.

def solve_chess_puzzle():
    """
    This function explains and prints the forced mate sequence.
    """
    print("White can force a checkmate in 5 moves.")
    print("The key line involves a queen check, followed by another check which sacrifices a knight but leads to an unstoppable mating attack.")
    print("Here is the winning sequence against optimal defense:")

    winning_sequence = [
        "1. Qxf7+ Kh8",
        "2. Qxe8+ Kh7",
        "3. Nf7+  Kg8",
        "4. Nh6++ Kh7", # This is a double check
        "5. Qg8#"
    ]

    for move in winning_sequence:
        print(move)

    num_white_moves = 5

    # Building the equation string as requested
    equation_parts = ["1"] * num_white_moves
    equation_string = " + ".join(equation_parts)

    print("\nEach of White's moves contributes one to the total number of moves to win.")
    # We are asked to output each number in the final equation.
    print(f"The calculation is: {equation_string} = {num_white_moves}")
    print(f"\nThus, White can win in {num_white_moves} moves.")

solve_chess_puzzle()
<<<5>>>