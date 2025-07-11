def solve_chess_puzzle():
    """
    This function prints the solution to the mate-in-2 chess puzzle.
    """
    move1_white = "Qd1"
    black_response = "fxf3"
    move2_white_mate = "Qd4#"
    
    print("The solution for the mate in 2 is:")
    print(f"1. {move1_white} {black_response}")
    print(f"2. {move2_white_mate}")

solve_chess_puzzle()