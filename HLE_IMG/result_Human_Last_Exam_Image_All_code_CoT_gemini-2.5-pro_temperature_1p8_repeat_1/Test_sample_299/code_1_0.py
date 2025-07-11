def solve_chess_puzzle():
    """
    This function explains and prints the solution to the provided mate-in-2 chess puzzle.
    """

    # White's first move (The key move)
    # The algebraic notation is Rf5. From square f3 to f5.
    white_move_1 = {"piece": "Rook", "from_rank": "3", "from_file": "f", "to_rank": "5", "to_file": "f"}

    # Black's response (main variation)
    # The algebraic notation is Qxf5.
    black_move_1 = {"piece": "Queen", "action": "takes", "target_square": "f5"}
    
    # White's second move (The mating move)
    # The algebraic notation is Ng6#. From square e7 to g6.
    white_move_2 = {"piece": "Knight", "from_rank": "7", "from_file": "e", "to_rank": "6", "to_file": "g"}

    print("The solution for the mate in 2 is as follows:")
    
    # Printing the first move
    print("\n1. White's first move is {} from {}{} to {}{}".format(
        white_move_1["piece"],
        white_move_1["from_file"],
        white_move_1["from_rank"],
        white_move_1["to_file"],
        white_move_1["to_rank"]
    ))

    # Printing Black's move
    print("   Black's main response is {} {}.".format(
        black_move_1["piece"],
        black_move_1["action"],
        black_move_1["target_square"]
    ))

    # Printing the final mating move with numbers from the coordinates separated as requested.
    print("\n2. White's second move delivers checkmate.")
    print("   The final move in the sequence is {} from {}{} to {}{}, checkmate.".format(
        white_move_2["piece"],
        white_move_2["from_file"],
        white_move_2["from_rank"],
        white_move_2["to_file"],
        white_move_2["to_rank"]
    ))
    
    print("\n--- Final Equation of moves ---")
    print("Move 1: R {} {} - {} {}".format(
        white_move_1["from_file"], 
        white_move_1["from_rank"], 
        white_move_1["to_file"], 
        white_move_1["to_rank"]
    ))
    print("Move 2: N {} {} - {} {} #".format(
        white_move_2["from_file"], 
        white_move_2["from_rank"], 
        white_move_2["to_file"], 
        white_move_2["to_rank"]
    ))

solve_chess_puzzle()