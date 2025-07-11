def solve_chess_puzzle():
    """
    This function explains and solves the chess puzzle.
    It prints the step-by-step mating sequence and the total number of moves.
    """
    
    # The puzzle is a mate-in-3 problem. Here is the move-by-move solution.
    
    move1_white = "Qxh7+"
    move1_black = "Rxh7"
    
    move2_white = "Ng6+"
    move2_black = "Kh8"
    
    move3_white = "Nf7#"
    
    # Step-by-step explanation
    print("White can deliver checkmate in 3 moves. Here is the sequence:")
    print("-" * 30)
    
    # Move 1
    print(f"1. {move1_white}")
    print("White sacrifices the Queen on h7. This is a check.")
    print("Black's King cannot capture due to the White Bishop on d3.")
    print(f"Black is forced to play: 1... {move1_black}\n")
    
    # Move 2
    print(f"2. {move2_white}")
    print("White moves the Knight from e4 to g6. This is a brilliant double check from both the Knight and the Bishop on d3.")
    print("In a double check, the King must move. The only legal square is h8.")
    print(f"Black is forced to play: 2... {move2_black}\n")
    
    # Move 3
    print(f"3. {move3_white}")
    print("White's Knight from e5 moves to f7, delivering checkmate.")
    print("The Black King on h8 is trapped. Its escape squares are controlled by the two White Knights.")
    print("-" * 30)
    
    # Final calculation
    moves_to_mate = 3
    print(f"The number of moves to mate the black king is {moves_to_mate}.")
    
    # Outputting each number in the final equation as requested
    print("The calculation is based on the number of moves for White:")
    print("1 (for Qxh7+) + 1 (for Ng6+) + 1 (for Nf7#) = 3")

solve_chess_puzzle()
<<<3>>>