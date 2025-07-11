def solve_chess_puzzle():
    """
    This function prints the solution to the mate-in-2 chess puzzle.
    """
    move1 = "1. Rxf4+"
    
    # First variation
    black_move1a = "1... Kxf4"
    white_mate1 = "2. Qd1#"
    
    # Second variation
    black_move1b = "1... Qf5"
    white_mate2 = "2. Re5#"
    
    print("The first move for White is a sacrifice:")
    print(move1)
    print("\nBlack has two possible replies. White delivers checkmate in either case.")
    print("\nVariation 1:")
    print(f"If Black plays {black_move1a}, White checkmates with {white_mate1}.")
    print("\nVariation 2:")
    print(f"If Black plays {black_move1b}, White checkmates with {white_mate2}.")

solve_chess_puzzle()