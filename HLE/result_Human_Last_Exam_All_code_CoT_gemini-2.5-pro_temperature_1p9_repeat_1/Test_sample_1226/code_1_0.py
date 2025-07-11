def solve_chess_puzzle():
    """
    Calculates the maximum material value for White in a mate-in-6 position.

    The position is constructed as follows:
    - A pawn at h2 will march to promote to a Queen, which takes 6 moves because
      a friendly Bishop at h4 must first move out of the way.
    - The Black King is caged at a1, with only b1 as a legal move.
    - All other squares, except for the path of the pawn and Bishop, are filled
      with the highest-value pieces (Queens) in a way that they cannot move
      or affect the outcome.
    """

    # Piece values
    queen_val = 9
    bishop_val = 3
    knight_val = 3
    pawn_val = 1

    # Count of each piece type based on the constructed mate-in-6 scenario
    
    # 47 filler Queens that are boxed in and cannot move
    num_queens = 47
    
    # 1 Bishop that must move out of the pawn's way on move 1
    num_bishops = 1
    
    # 1 Knight used to cage the Black King
    num_knights = 1
    
    # 1 Pawn for the promotion mate + 6 pawns to cage the Black King
    num_pawns = 1 + 6
    
    # Calculate the total points
    queen_points = num_queens * queen_val
    bishop_points = num_bishops * bishop_val
    knight_points = num_knights * knight_val
    pawn_points = num_pawns * pawn_val
    
    total_points = queen_points + bishop_points + knight_points + pawn_points

    # Output the final equation with each component
    print(f"The calculation for the total points is:")
    print(f"({num_queens} Queens * {queen_val} pts) + ({num_bishops} Bishop * {bishop_val} pts) + ({num_knights} Knight * {knight_val} pts) + ({num_pawns} Pawns * {pawn_val} pt)")
    print(f"= {queen_points} + {bishop_points} + {knight_points} + {pawn_points}")
    print(f"= {total_points}")
    print("\nThe greatest number of points is:")
    print(total_points)

solve_chess_puzzle()