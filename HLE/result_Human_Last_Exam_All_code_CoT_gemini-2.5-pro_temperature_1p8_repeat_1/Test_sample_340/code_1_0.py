import sys
# This script explains the solution to the Three-check chess puzzle.

# The puzzle asks for the minimal number of moves for White to win from the given position.
# A win is achieved by delivering 3 checks or by checkmate.
# The initial FEN indicates the check count is 0 for both sides.

# We will demonstrate the shortest forced winning line for White.

def solve_three_check_puzzle():
    """
    This function prints the step-by-step solution to the chess puzzle.
    """
    print("The minimal number of moves for White to win is 3.")
    print("The winning strategy involves delivering three consecutive checks.\n")
    
    # Define the winning sequence of moves
    winning_moves = [
        {"move": "1. Qe6+", "response": "1... Be7", "check_number": 1},
        {"move": "2. Qxe7+", "response": "2... Kxe7", "check_number": 2},
        {"move": "3. O-O-O+", "response": "White wins", "check_number": 3},
    ]
    
    white_move_count = 0
    
    for item in winning_moves:
        white_move_count += 1
        print(f"White's Move {white_move_count}: {item['move']}")
        print(f"This delivers check number {item['check_number']} to the Black king.")
        if item['response'] == "White wins":
            print("With the third check delivered, White wins the game.")
        else:
            print(f"Black's forced response is {item['response']}.")
        print("-" * 30)

    print(f"\nFinal conclusion: The minimal number of moves for White to win is {white_move_count}.")

if __name__ == '__main__':
    solve_three_check_puzzle()
