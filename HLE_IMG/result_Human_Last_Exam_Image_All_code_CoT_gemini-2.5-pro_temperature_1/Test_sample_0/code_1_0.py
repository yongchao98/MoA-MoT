def solve_chess_puzzle():
    """
    This function prints the solution to the chess puzzle.
    The puzzle is a mate in 2 for Black, without moving the queens.
    """
    move1 = "Ne3+"
    move2 = "Nf1#"
    
    print("The mate-in-2 sequence for Black is:")
    print(f"1. {move1}")
    print(f"2. {move2}")

solve_chess_puzzle()