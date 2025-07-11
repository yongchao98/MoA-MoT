import collections

def solve_chess_puzzle():
    """
    This function explains and prints the solution to the chess puzzle.
    """
    print("The solution involves a classic sacrificial attack against the Black king.")
    print("White can force a checkmate in 4 moves against Black's best defense.")
    print("\nHere is the winning sequence:\n")

    # Using a list of tuples to store moves and their explanations
    # (White's move, Black's move, Explanation)
    moves = [
        ("Nxg6", "hxg6", "White sacrifices the knight on e5 to destroy the pawn cover around the Black king. Black's best defense is to capture the knight, as other moves like ...Qe8 or ...Kh8 lead to a faster mate (in 2 or 3 moves respectively)."),
        ("Qxg6+", "Kh8", "White's Queen takes the now-undefended g6 pawn, delivering a check. The Black king has only one legal move, retreating to the corner."),
        ("Ng5", "Qf8", "White brings the second knight from e4 into the attack. This creates an unstoppable checkmate threat on h7. Black has no move to prevent this. We'll use ...Qf8 as an example of a move that doesn't change the outcome."),
        ("Qh7#", "", "White's Queen delivers checkmate on h7, protected by the knight on g5. The Black king has no escape squares.")
    ]

    for i, (white_move, black_move, explanation) in enumerate(moves, 1):
        move_number = str(i)
        
        # Format the output clearly
        if black_move:
            print(f"{move_number}. {white_move} {black_move}")
        else:
            print(f"{move_number}. {white_move}")
            
        print(f"   - {explanation}\n")
    
    final_answer = len(moves)
    print(f"Conclusion: It takes White {final_answer} moves to deliver checkmate.")

solve_chess_puzzle()