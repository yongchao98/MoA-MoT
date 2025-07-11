def solve_three_check_puzzle():
    """
    This function explains the solution to the three-check chess puzzle.
    It prints the minimal number of moves for White to achieve a forced win.
    """
    
    # The line of play demonstrates the forced win, assuming optimal defense from Black.
    # The length of a forced win is determined by the opponent's best defense (the one that prolongs the game the most).
    
    line_of_play = [
        "1. Bxd7+ (Check 1)  Qxd7",
        "2. Qxe5+ (Check 2)  Be7",
        "3. O-O-O           O-O (Black's best defense, prolonging the game)",
        "4. Rxd7            Qxd7 (Black's best defense in the sub-variation)",
        "5. Rxd7            Nxd7",
        "6. Qxe7+ (Check 3)  ...",
    ]
    
    winning_moves_by_white = 6
    
    print("The shortest forced winning line for White is as follows:")
    for move in line_of_play:
        print(move)
        
    print("\nWhite delivers the third check on the 6th move.")
    print(f"Minimal amount of moves by white to win: {winning_moves_by_white}")

solve_three_check_puzzle()