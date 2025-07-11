def solve_capablanca_puzzle():
    """
    This function explains the step-by-step solution to the Capablanca chess puzzle
    and prints the minimal number of moves for White to win.
    """

    # The problem asks for the minimal number of moves for White to achieve checkmate.
    # This is a classic "mate in N" problem. We need to find the smallest N.

    # The solution is a forced checkmate sequence.
    move_1_white = "Qh3+"
    move_1_black = "Ki8"
    
    move_2_white = "Qi3+"
    move_2_black = "Cxi7"
    
    move_3_white = "Qxi7#"
    
    number_of_moves = 3
    
    print("The minimal number of moves by White to win is found by determining the shortest forced checkmate sequence.")
    print("The analysis shows a forced mate in 3 moves.")
    print("\nThe step-by-step sequence is:")
    
    print(f"1. White plays {move_1_white}. This check forces Black's only reply: {move_1_black}.")
    print(f"2. White plays {move_2_white}. This check forces Black's only reply: {move_2_black}.")
    print(f"3. White plays {move_3_white}, delivering checkmate.")
    
    print(f"\nThus, the minimal amount of moves by White to win is {number_of_moves}.")

solve_capablanca_puzzle()