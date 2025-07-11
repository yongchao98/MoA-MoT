def calculate_max_material():
    """
    Calculates the maximum material points for a mate in 6 moves.

    This is based on a "self-block" construction where White must spend
    five moves clearing a path for a mating piece.

    1.  A White King is on the board.
    2.  The Black King needs 2 squares to shuttle between (e.g., a1, b1).
    3.  A mating attack requires a path to be cleared. In our construction, this
        means 5 queens must move out of the way, and the mating piece needs
        a path. This requires 7 additional empty squares.
        - 5 squares for the blocking pieces to move to.
        - 2 squares for the mating piece's final path.
    4.  The total squares unavailable for placing material are:
        - 1 for the White King.
        - 2 for the Black King's shuttle area.
        - 7 for the empty squares needed by the mating mechanism.
    5.  The remaining squares can be filled with Queens for maximum points.
    """
    total_squares = 64
    queen_points = 9

    # Squares that cannot be filled with scoring material
    white_king_squares = 1
    black_king_shuttle_area = 2 # a1, b1
    empty_squares_for_mechanism = 7 # 5 for queen moves + 2 for mating path

    unavailable_squares = white_king_squares + black_king_shuttle_area + empty_squares_for_mechanism

    # The rest of the board can be filled with Queens.
    # This count includes the queens that are part of the cage/mating mechanism.
    number_of_queens = total_squares - unavailable_squares

    total_points = number_of_queens * queen_points

    print("This problem is a specific example of a chess composition theme called a 'maximummer'.")
    print("The solution involves a 'self-block' mechanism where White must spend 5 moves clearing a path for a mate on the 6th move.")
    print(f"To achieve this, certain squares must be left empty for the mechanism to work, and one square is occupied by the White King.")
    print(f"Number of squares on a chessboard: {total_squares}")
    print(f"Squares unavailable for material (White King + Black King Area + Empty Mechanism Squares): {unavailable_squares}")
    print(f"Number of squares available for White Queens: {number_of_queens}")
    print("\nFinal point calculation:")
    print(f"{number_of_queens} (Queens) * {queen_points} (points per Queen) = {total_points}")

calculate_max_material()