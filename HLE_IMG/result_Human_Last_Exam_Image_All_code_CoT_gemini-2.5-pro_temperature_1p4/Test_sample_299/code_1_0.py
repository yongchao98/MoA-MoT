def solve_chess_puzzle():
    """
    This script prints the solution to the chess mate-in-two puzzle.
    """
    move1 = "1. Qf4+"
    
    # First variation
    black_reply1 = "1... exf4"
    white_mate1 = "2. Nd6#"
    
    # Second variation
    black_reply2 = "1... Kd3"
    white_mate2 = "2. Qe4#"
    
    print("The solution to the mate-in-two puzzle is as follows:")
    print(f"White's first move: {move1}")
    print("\nThere are two possible responses from Black:")
    
    print("\nVariation 1:")
    print(f"If Black plays {black_reply1}, White checkmates with {white_mate1}.")
    print("1. Q a4 -> f4 +")
    print("1. ... p e5 x f4")
    print("2. N e7 -> d6 #")


    print("\nVariation 2:")
    print(f"If Black plays {black_reply2}, White checkmates with {white_mate2}.")
    print("1. Q a4 -> f4 +")
    print("1. ... K e4 -> d3")
    print("2. Q f4 -> e4 #")

solve_chess_puzzle()
