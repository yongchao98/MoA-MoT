def solve_chess_puzzle():
    """
    This function prints the solution to the chess puzzle.
    The puzzle involves a forced checkmate for Black.
    The solution is a classic smothered mate sequence.
    """
    print("The best sequence of moves that forces checkmate is:")
    
    # The sequence of moves for the smothered mate
    move_1_black = "1... Qg1+"
    move_2_white = "2. Rxg1"
    move_2_black = "2... Nf2#" # This is the second move for black, which happens on the second turn.

    print(move_1_black)
    print(move_2_white)
    print(move_2_black)

solve_chess_puzzle()