def solve_chess_puzzle():
    """
    This function explains and prints the solution to the chess puzzle.
    The puzzle is a mate in 3 for White.
    """
    print("The solution is a forced checkmate in 3 moves for White.")
    print("Here is the winning sequence and the explanation for the main line:")
    print("-" * 30)

    # Define the moves of the main mating line
    move1_white = "Nxf7+"
    move1_black = "Rxf7"
    move2_white = "Ng6+"
    move2_black = "hxg6"
    move3_white = "Qh8#"

    # Print the "equation" of moves
    print(f"1. {move1_white} {move1_black}")
    print("Explanation: White's Knight captures the pawn, delivering check. Black's only legal move is to recapture with the Rook.")
    print("")

    print(f"2. {move2_white} {move2_black}")
    print("Explanation: White's other Knight moves to create a discovered check from the Queen. Black is forced to capture the Knight.")
    print("")

    print(f"3. {move3_white}")
    print("Explanation: White's Queen delivers the final checkmate.")
    print("-" * 30)
    
    number_of_moves = 3
    print(f"Total moves to mate: {number_of_moves}")

solve_chess_puzzle()