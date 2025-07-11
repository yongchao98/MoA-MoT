def solve_chess_puzzle():
    """
    This function prints the best sequence of moves for black to checkmate white.
    The sequence is known as a "smothered mate".
    """
    # 1... Qg1+: Black sacrifices the queen, forcing the rook to capture.
    # 2. Rxg1 Nf2#: Black's knight delivers the final checkmate blow.
    
    print("The best sequence of moves is:")
    print("1... Qg1+")
    print("2. Rxg1 Nf2#")

solve_chess_puzzle()