import re

def solve_chess_puzzle():
    """
    Analyzes the given chess position and provides the best move for White.
    """
    # FEN string for the position.
    fen = "8/3p4/1kpP4/p1q5/P7/8/5Q2/6K1 w - - 0 1"
    
    # The best move for White, determined by strategic analysis.
    best_move = "Qxc5+"
    
    explanation = (
        "The optimal move for White is Qxc5+. This forces a queen trade, "
        "transitioning the game into a pawn endgame. In this endgame, White's "
        "passed pawn on d6, combined with the ability for the White king to activate, "
        "creates a decisive and winning advantage for White."
    )

    print(f"Position FEN: {fen}")
    print(f"White to move. Best move analysis:")
    print(explanation)
    
    print("\nThe best move is:")
    print(best_move)
    
    # The prompt requires outputting each number in the 'final equation' (the move notation).
    # The following code finds and prints all digits present in the move string.
    numbers_in_move = re.findall(r'\d', best_move)
    
    print("\nPrinting each number found in the move notation:")
    if numbers_in_move:
        for number in numbers_in_move:
            print(number)
    else:
        # This case would handle moves like 'e4' which have no numbers.
        print("No numbers were found in this move's notation.")

# Execute the function to solve the puzzle.
solve_chess_puzzle()

<<<Qxc5+>>>