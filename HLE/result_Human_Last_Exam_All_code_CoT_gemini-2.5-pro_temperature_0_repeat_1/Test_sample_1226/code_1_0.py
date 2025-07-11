def calculate_max_chess_material():
    """
    This script calculates the greatest number of points of white material
    on a chessboard where white can mate, but not in less than 6 moves.
    """

    # Define the constants for the problem
    total_board_squares = 64
    points_per_queen = 9
    number_of_kings = 2
    minimum_mate_duration = 6

    # Explain the logic
    print("Step 1: To maximize the material score, we must use the most valuable piece, the Queen (9 points).")
    print("The goal is to fit as many queens as possible onto the 64 squares of the board.")
    print("-" * 20)

    print("Step 2: Account for required pieces and empty squares.")
    print(f"A White King and a Black King must be on the board, occupying {number_of_kings} squares.")
    
    # Determine the number of empty squares needed to delay the mate.
    # A mate in M moves requires the black king to make M-1 moves into empty squares.
    # To meet the "not less than 6 moves" criteria and maximize queens, we aim for a mate in exactly 6 moves.
    required_king_moves = minimum_mate_duration - 1
    print(f"To delay the mate to {minimum_mate_duration} moves, the Black King needs an escape path of {required_king_moves} empty squares.")
    print("-" * 20)

    # Calculate the number of queens
    print("Step 3: Calculate the maximum number of queens.")
    non_queen_squares = number_of_kings + required_king_moves
    print(f"The total number of squares that cannot be queens is {number_of_kings} (for kings) + {required_king_moves} (empty) = {non_queen_squares}.")
    
    number_of_queens = total_board_squares - non_queen_squares
    print(f"The maximum number of queens is {total_board_squares} (total) - {non_queen_squares} (non-queen) = {number_of_queens}.")
    print("-" * 20)

    # Calculate the final score
    print("Step 4: Calculate the final point value.")
    total_points = number_of_queens * points_per_queen
    
    print("Final Equation:")
    print(f"{number_of_queens} Queens * {points_per_queen} points/Queen = {total_points}")


calculate_max_chess_material()
<<<513>>>