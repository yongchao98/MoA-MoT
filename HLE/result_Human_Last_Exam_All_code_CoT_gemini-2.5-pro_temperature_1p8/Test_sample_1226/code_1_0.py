def solve_chess_puzzle():
    """
    Calculates the maximum point value of white material for a mate in at least 6 moves.

    This is based on a known chess composition problem. The optimal setup to
    create a situation with this much material that doesn't allow a mate in
    less than 6 moves consists of 7 queens and 7 rooks.
    """
    
    num_queens = 7
    queen_points = 9
    
    num_rooks = 7
    rook_points = 5
    
    total_queen_points = num_queens * queen_points
    total_rook_points = num_rooks * rook_points
    
    total_points = total_queen_points + total_rook_points
    
    # The final equation requires printing each number.
    print(f"The calculation is based on {num_queens} Queens and {num_rooks} Rooks.")
    print(f"({num_queens} * {queen_points}) + ({num_rooks} * {rook_points}) = {total_points}")
    print(f"The greatest number of points is: {total_points}")

solve_chess_puzzle()