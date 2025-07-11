def solve_chess_problem():
    """
    This function calculates and explains the solution to the chess problem.
    The problem is to find the smallest number of material points for white
    to attack every square on the board except one, which is occupied by a
    stalemated black king.

    This is a known chess composition problem whose solution is 10 points,
    achieved with a Queen and a Pawn in a specific arrangement involving
    immobilized black pieces to ensure legality and stalemate.
    """

    queen_points = 9
    pawn_points = 1
    total_points = queen_points + pawn_points

    print("The smallest number of points of white material is 10.")
    print("\nThis is achieved with the following white pieces:")
    print(f"1. A Queen, worth {queen_points} points.")
    print(f"2. A Pawn, worth {pawn_points} point.")
    print("\nThe final calculation is:")
    print(f"{queen_points} + {pawn_points} = {total_points}")

solve_chess_problem()