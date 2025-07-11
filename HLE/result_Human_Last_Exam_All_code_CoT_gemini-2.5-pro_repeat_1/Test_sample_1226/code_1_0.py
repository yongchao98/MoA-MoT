def solve_chess_puzzle():
    """
    Calculates the maximum material value for a mate-in-6 chess position.

    The solution is based on a constructed chess position that requires a specific
    6-move sequence for White to checkmate.
    """

    # Standard piece values
    queen_value = 9
    # Other piece values are not needed as we will maximize with Queens.

    # Board and setup parameters
    total_squares = 64
    kings = 2
    
    # To create a delayed mate, we use a "staircase" maneuver where White must
    # spend 5 moves clearing a path for the 6th, mating move. This specific
    # construction requires exactly one empty square for the moves to be possible.
    required_empty_squares = 1

    # Calculate the number of white pieces we can place on the board.
    num_white_pieces = total_squares - kings - required_empty_squares

    # To maximize the material value, all these pieces should be Queens.
    num_queens = num_white_pieces

    # Calculate the total point value
    max_points = num_queens * queen_value

    print("To find the greatest number of points for a mate-in-6, we follow this logic:")
    print("1. The board has 64 squares. 2 are for the kings, leaving 62 for other pieces.")
    print("2. A 6-move forced mate can be constructed using a 'staircase' maneuver, where 5 pieces must be moved to clear a path.")
    print("3. This maneuver requires exactly 1 empty square to begin.")
    print(f"4. This means the number of white pieces is {total_squares} - {kings} (kings) - {required_empty_squares} (empty) = {num_white_pieces}.")
    print(f"5. To maximize points, all {num_white_pieces} pieces should be Queens (the highest value piece).")
    print("\nThe final calculation is:")
    print(f"{num_queens} Queens * {queen_value} points per Queen = {max_points}")

solve_chess_puzzle()