def solve_chess_puzzle():
    """
    Calculates the maximum point value of white material for a mate in >= 6 moves.

    The solution is based on constructing a position with two parts:
    1. A large "paralyzed" grid of high-value pieces that cannot move.
    2. A smaller area where a slow checkmating combination (mate in >= 6) occurs.

    - Paralyzed Grid (files a-f):
      - We use alternating rows of Queens and Rooks.
      - 4 rows of Queens * 6 files = 24 Queens.
      - 4 rows of Rooks * 6 files = 24 Rooks.
      - These pieces are all immobilized by each other.

    - Mating Force (files g-h):
      - A force of 3 Rooks is chosen. It's powerful enough to guarantee a mate against a lone king
        but slow enough that a mate-in-6-or-more position can be constructed.

    The total point value is the sum of all these pieces.
    """
    
    num_queens = 24
    queen_value = 9
    
    # Rooks are used in both the paralyzed grid and the mating force
    num_paralyzed_rooks = 24
    num_mating_rooks = 3
    rook_value = 5
    
    total_rooks = num_paralyzed_rooks + num_mating_rooks
    
    total_value = (num_queens * queen_value) + (total_rooks * rook_value)
    
    print(f"The calculation is based on the sum of points from Queens and Rooks.")
    print(f"There are {num_queens} Queens and {total_rooks} Rooks.")
    print(f"Point values: Queen=9, Rook=5.")
    print(f"Equation: ({num_queens} * {queen_value}) + ({total_rooks} * {rook_value}) = {total_value}")
    print(f"Greatest number of points: {total_value}")

solve_chess_puzzle()
