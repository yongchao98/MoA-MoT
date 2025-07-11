def solve_chess_puzzle():
    """
    This function prints the step-by-step moves for the forced checkmate.
    """
    moves = [
        "1... Qg1+",
        "2. Rxg1",
        "2... Nf2#"
    ]
    print("The best sequence of moves that forces a checkmate is:")
    for move in moves:
        print(move)

solve_chess_puzzle()