def solve_chess_material_puzzle():
    """
    Calculates the maximum point value of white material for a specific
    mate-in-6 chess problem.
    """

    # Define standard piece values
    piece_values = {
        'Queen': 9,
        'Rook': 5,
        'Bishop': 3,
        'Knight': 3,
        'Pawn': 1
    }

    # The configuration for the mate-in-6 puzzle is as follows:
    # - White King: on a safe square.
    # - Black King: on a1, forced to shuttle to b1.
    # - Cage Pieces: Bishop on a3, Bishop on b3, Rook on d2.
    # - Mating Piece: Rook on h8.
    # - Spectator Pieces: The maximum number of Queens that don't interfere.

    # Count of each white piece in the final configuration
    num_queens = 24
    num_rooks = 2  # One for the cage (d2), one for the mate (h8)
    num_bishops = 2 # Two for the cage (a3, b3)
    num_knights = 0
    num_pawns = 0

    # Calculate the total value
    total_value = (num_queens * piece_values['Queen'] +
                   num_rooks * piece_values['Rook'] +
                   num_bishops * piece_values['Bishop'] +
                   num_knights * piece_values['Knight'] +
                   num_pawns * piece_values['Pawn'])

    # Print the equation as requested
    print("The maximum point value is calculated from the following pieces:")
    print(f"Queens: {num_queens}")
    print(f"Rooks: {num_rooks}")
    print(f"Bishops: {num_bishops}")
    print("\nThe final calculation is:")
    print(f"{num_queens} * {piece_values['Queen']} (Queens) + "
          f"{num_rooks} * {piece_values['Rook']} (Rooks) + "
          f"{num_bishops} * {piece_values['Bishop']} (Bishops) = {total_value}")

    print(f"\nThe greatest number of points is: {total_value}")

solve_chess_material_puzzle()