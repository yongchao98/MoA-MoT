def solve_chess_mate_in_2():
    """
    This function calculates and prints the solution to the given chess puzzle.
    The solution is a mate in 2 for Black, without moving the queens.
    The code explicitly uses the numbers from the chess notation to form the output string.
    """
    
    # The first move is Rook to f2, delivering a check.
    # The number in this move's coordinate is 2.
    move1_rank = 2
    
    # The second move is Bishop to f1, delivering checkmate.
    # The number in this move's coordinate is 1.
    move2_rank = 1

    # Construct the standard algebraic notation for the moves.
    first_move = f"Rf{move1_rank}+"
    second_move = f"Bf{move2_rank}#"
    
    # Print the final sequence of moves for Black.
    print(f"The mate-in-2 sequence for Black is: {first_move}, {second_move}")

solve_chess_mate_in_2()